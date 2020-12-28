use crate::simplification_common::*;
use crate::tables::*;
use crate::ForrusttsError;
use crate::SamplesInfo;
use crate::SimplificationBuffers;
use crate::SimplificationFlags;
use crate::SimplificationOutput;

/// Simplify a [``TableCollection``].
///
/// # Parameters
///
/// * `samples`:
/// * `flags`: modify the behavior of the simplification algorithm.
/// * `tables`: a [``TableCollection``] to simplify.
/// * `output`: Where simplification output gets written.
///             See [``SimplificationOutput``].
///
/// # Notes
///
/// The input tables must be sorted.
/// See [``TableCollection::sort_tables_for_simplification``].
///
/// It is common to simplify many times during a simulation.
/// To avoid making big allocations each time, see
/// [``simplify_tables``] to keep memory allocations
/// persistent between simplifications.
pub fn simplify_tables_without_state(
    samples: &SamplesInfo,
    flags: SimplificationFlags,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) -> Result<(), ForrusttsError> {
    let mut state = SimplificationBuffers::new();
    simplify_tables(samples, flags, &mut state, tables, output)
}

/// Simplify a [``TableCollection``].
///
/// This differs from [``simplify_tables_without_state``] in that the big memory
/// allocations made during simplification are preserved in
/// an instance of [``SimplificationBuffers``].
///
/// # Parameters
///
/// * `samples`:
/// * `flags`: modify the behavior of the simplification algorithm.
/// * `state`: These are the internal data structures used
///            by the simpilfication algorithm.
/// * `tables`: a [``TableCollection``] to simplify.
/// * `output`: Where simplification output gets written.
///             See [``SimplificationOutput``].
///
/// # Notes
///
/// The input tables must be sorted.
/// See [``TableCollection::sort_tables_for_simplification``].
pub fn simplify_tables(
    samples: &SamplesInfo,
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

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: we need lots more tests of these validations!

    #[test]
    fn test_simplify_tables_unsorted_edges() {
        let mut tables = TableCollection::new(1000).unwrap();

        tables.add_node(0, 0).unwrap(); // parent
        tables.add_node(1, 0).unwrap(); // child
        tables.add_edge(100, tables.genome_length(), 0, 1).unwrap();
        tables.add_edge(0, 100, 0, 1).unwrap();

        let mut output = SimplificationOutput::new();

        let mut samples = SamplesInfo::new();
        samples.samples.push(1);

        let _ = simplify_tables_without_state(
            &samples,
            SimplificationFlags::VALIDATE_ALL,
            &mut tables,
            &mut output,
        )
        .map_or_else(
            |x: ForrusttsError| {
                assert_eq!(
                    x,
                    ForrusttsError::TablesError {
                        value: TablesError::EdgesNotSortedByLeft
                    }
                )
            },
            |_| panic!(),
        );
    }
}
