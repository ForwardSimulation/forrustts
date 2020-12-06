use crate::nested_forward_list::NestedForwardList;
use crate::tables::*;
use crate::tsdef::{SamplesVec, TsInt, NULL};

fn swap_edges(tables: &mut TableCollection, edges: &mut EdgeTable) {
    std::mem::swap(&mut tables.edges_, edges);
}

fn swap_nodes(tables: &mut TableCollection, nodes: &mut NodeTable) {
    std::mem::swap(&mut tables.nodes_, nodes);
}

fn setup_idmap(nodes: &NodeTable) -> SamplesVec {
    return vec![NULL; nodes.len()];
}

struct Segment {
    pub left: i64,
    pub right: i64,
    pub node: TsInt,
}

type AncestryList = NestedForwardList<Segment>;

pub fn simplify_tables(samples: &SamplesVec, tables: &mut TableCollection) -> SamplesVec {
    if tables.sites_.len() > 0 || tables.mutations_.len() > 0 {
        panic!("mutation simplification not yet implemented");
    }

    let idmap = setup_idmap(&tables.nodes_);
    let mut new_nodes = NodeTable::new();
    let mut new_edges = EdgeTable::new();
    let mut ancestry = AncestryList::new();

    swap_edges(tables, &mut new_edges);
    swap_nodes(tables, &mut new_nodes);
    return idmap;
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_swap_edges() {
        let mut tables = create_table_collection(1000);
        let mut edges = EdgeTable::new();

        let num_edges = edge_table_add_row(&mut edges, 0, 1, 3, 4).unwrap();
        assert_eq!(1, num_edges);

        swap_edges(&mut tables, &mut edges);

        assert_eq!(0, edges.len());
        assert_eq!(num_edges, tables.edges_.len());
    }
}
