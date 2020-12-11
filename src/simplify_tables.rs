use crate::simplification_buffers::SimplificationBuffers;
use crate::simplification_logic;
use crate::tables::*;
use crate::tsdef::SamplesVec;

pub fn simplify_tables(samples: &SamplesVec, tables: &mut TableCollection) -> SamplesVec {
    let mut state = SimplificationBuffers::new();
    return simplify_tables_with_state(samples, &mut state, tables);
}

pub fn simplify_tables_with_state(
    samples: &SamplesVec,
    state: &mut SimplificationBuffers,
    tables: &mut TableCollection,
) -> SamplesVec {
    if tables.sites_.len() > 0 || tables.mutations_.len() > 0 {
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
            &mut state.temp_edge_buffer,
            &mut state.new_nodes,
            &mut state.new_edges,
            &mut state.ancestry,
            &mut state.overlapper,
            &mut idmap,
        );

        if state.new_edges.len() >= 1024 && new_edges_inserted + state.new_edges.len() < edge_i {
            for i in 0..state.new_edges.len() {
                tables.edges_[new_edges_inserted + i] = state.new_edges[i];
            }
            new_edges_inserted += state.new_edges.len();
            state.new_edges.clear();
        }
    }

    if state.new_edges.len() > 0 {
        for i in 0..state.new_edges.len() {
            tables.edges_[new_edges_inserted + i] = state.new_edges[i];
        }
        new_edges_inserted += state.new_edges.len();
        state.new_edges.clear();
    }

    tables.edges_.drain(new_edges_inserted..tables.edges_.len());
    std::mem::swap(&mut tables.nodes_, &mut state.new_nodes);

    return idmap;
}
