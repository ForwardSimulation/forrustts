use crate::simplification_logic;
use crate::tables::*;
use crate::tsdef::SamplesVec;

pub fn simplify_tables(samples: &SamplesVec, tables: &mut TableCollection) -> SamplesVec {
    if tables.sites_.len() > 0 || tables.mutations_.len() > 0 {
        panic!("mutation simplification not yet implemented");
    }

    let mut idmap = simplification_logic::setup_idmap(&tables.nodes_);
    let mut new_nodes = NodeTable::new();
    let mut temp_edge_buffer = EdgeTable::new();
    let mut new_edges = EdgeTable::new();
    let mut ancestry = simplification_logic::AncestryList::new();
    let mut overlapper = simplification_logic::SegmentOverlapper::new();

    ancestry.reset(tables.num_nodes());

    simplification_logic::record_sample_nodes(
        &samples,
        &tables,
        &mut new_nodes,
        &mut ancestry,
        &mut idmap,
    );

    let mut edge_i = 0;
    let num_edges = tables.num_edges();
    while edge_i < num_edges {
        let u = tables.edges_[edge_i].parent;
        edge_i = simplification_logic::find_parent_child_segment_overlap(
            &tables.edges_,
            edge_i,
            num_edges,
            tables.get_length(),
            u,
            &mut ancestry,
            &mut overlapper,
        );

        simplification_logic::merge_ancestors(
            &tables.nodes_,
            tables.get_length(),
            u,
            &mut temp_edge_buffer,
            &mut new_nodes,
            &mut new_edges,
            &mut ancestry,
            &mut overlapper,
            &mut idmap,
        );
    }
    simplification_logic::swap_edges(tables, &mut new_edges);
    simplification_logic::swap_nodes(tables, &mut new_nodes);

    return idmap;
}
