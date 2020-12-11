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
    let mut new_edges_inserted: usize = 0;
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

        if new_edges.len() >= 1024 && new_edges_inserted + new_edges.len() < edge_i {
            for i in 0..new_edges.len() {
                tables.edges_[new_edges_inserted + i] = new_edges[i];
            }
            new_edges_inserted += new_edges.len();
            new_edges.clear();
        }
    }

    if new_edges.len() > 0 {
        for i in 0..new_edges.len() {
            tables.edges_[new_edges_inserted + i] = new_edges[i];
        }
        new_edges_inserted += new_edges.len();
        new_edges.clear();
    }

    tables.edges_.drain(new_edges_inserted..tables.edges_.len());
    simplification_logic::swap_nodes(tables, &mut new_nodes);

    return idmap;
}
