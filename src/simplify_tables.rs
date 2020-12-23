use crate::simplification_common::*;
use crate::tables::*;
use crate::ForrusttsError;
use crate::IdType;
use crate::SimplificationBuffers;
use crate::SimplificationFlags;
use crate::SimplificationOutput;

pub fn simplify_tables(
    samples: &[IdType],
    flags: SimplificationFlags,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) -> Result<(), ForrusttsError> {
    let mut state = SimplificationBuffers::new();
    simplify_tables_with_buffers(samples, flags, &mut state, tables, output)
}

pub fn simplify_tables_with_buffers(
    samples: &[IdType],
    flags: SimplificationFlags,
    state: &mut SimplificationBuffers,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) -> Result<(), ForrusttsError> {
    setup_simplification(samples, tables, flags, state, output)?;

    let mut edge_i = 0;
    let num_edges = tables.num_edges();
    let mut new_edges_inserted: usize = 0;
    while edge_i < num_edges {
        edge_i = process_parent(
            tables.edges_[edge_i].parent,
            (edge_i, num_edges),
            &tables,
            state,
            output,
        )?;

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

    Ok(())
}
