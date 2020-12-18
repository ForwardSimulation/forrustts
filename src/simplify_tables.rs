use crate::simplification_buffers::SimplificationBuffers;
use crate::simplification_logic;
use crate::tables::*;
use crate::IdType;

pub fn simplify_tables(samples: &[IdType], tables: &mut TableCollection) -> Vec<IdType> {
    let mut state = SimplificationBuffers::new();
    simplify_tables_with_buffers(samples, &mut state, tables)
}

pub fn simplify_tables_with_buffers(
    samples: &[IdType],
    state: &mut SimplificationBuffers,
    tables: &mut TableCollection,
) -> Vec<IdType> {
    if !tables.sites_.is_empty() || !tables.mutations_.is_empty() {
        panic!("mutation simplification not yet implemented");
    }

    let mut idmap = simplification_logic::setup_idmap(&tables.nodes_);

    state.clear();
    state.ancestry.reset(tables.num_nodes());

    simplification_logic::record_sample_nodes(
        &samples,
        &tables,
        &mut state.new_nodes,
        &mut state.ancestry,
        &mut idmap,
    );

    let mut edge_i = 0;
    let num_edges = tables.num_edges();
    let mut new_edges_inserted: usize = 0;
    while edge_i < num_edges {
        let u = tables.edges_[edge_i].parent;
        edge_i = simplification_logic::find_parent_child_segment_overlap(
            &tables.edges_,
            edge_i,
            num_edges,
            tables.get_length(),
            u,
            &mut state.ancestry,
            &mut state.overlapper,
        );

        simplification_logic::merge_ancestors(
            &tables.nodes_,
            tables.get_length(),
            u,
            state,
            &mut idmap,
        );

        if state.new_edges.len() >= 1024 && new_edges_inserted + state.new_edges.len() < edge_i {
            for i in state.new_edges.drain(..) {
                tables.edges_[new_edges_inserted] = i;
                new_edges_inserted += 1;
            }
            assert_eq!(state.new_edges.len(), 0);
        }
    }

    tables.edges_.truncate(new_edges_inserted);
    tables.edges_.append(&mut state.new_edges);
    std::mem::swap(&mut tables.nodes_, &mut state.new_nodes);

    idmap
}
