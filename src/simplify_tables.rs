use crate::simplification_logic;
use crate::tables::*;
use crate::IdType;
use crate::SimplificationBuffers;
use crate::SimplificationFlags;
use crate::SimplificationOutput;

pub fn simplify_tables(
    samples: &[IdType],
    flags: SimplificationFlags,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) {
    let mut state = SimplificationBuffers::new();
    simplify_tables_with_buffers(samples, flags, &mut state, tables, output)
}

pub fn simplify_tables_with_buffers(
    samples: &[IdType],
    flags: SimplificationFlags,
    state: &mut SimplificationBuffers,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) {
    if !tables.sites_.is_empty() || !tables.mutations_.is_empty() {
        panic!("mutation simplification not yet implemented");
    }

    if flags.bits() != 0 {
        panic!("SimplificationFlags must be zero");
    }

    simplification_logic::setup_idmap(&tables.nodes_, &mut output.idmap);

    state.clear();
    state.ancestry.reset(tables.num_nodes());

    simplification_logic::record_sample_nodes(
        &samples,
        &tables,
        &mut state.new_nodes,
        &mut state.ancestry,
        &mut output.idmap,
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
            &mut output.idmap,
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
}
