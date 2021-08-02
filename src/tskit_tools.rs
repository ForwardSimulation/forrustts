//! Data interchange to ``tskit`` format using [``tskit``].
//!
//! # Note
//!
//! This module may be prone to API-breaking changes!
//! Currently, the function implementations are the minimum
//! needed for testing/data exchange.
//!
//! As things develop, we may either add new functions or
//! refactor the existing.

use crate::newtypes::Time;
use crate::TableCollection;
use tskit::{tsk_flags_t, tsk_id_t, TSK_NODE_IS_SAMPLE, TSK_NULL};

/// Return a closure to help reverse time.
///
/// For all input values, ``t`` the closure will
/// return ``-1.0*(t - x) as f64``.
pub fn simple_time_reverser(x: Time) -> Box<dyn Fn(Time) -> f64> {
    Box::new(move |t: Time| -1. * (t.0 - x.0) as f64)
}

/// Convert a [``TableCollection``](crate::TableCollection)
/// to ``tskit`` format.
///
/// # Parameters
///
/// * `tables`: A [``TableCollection``](crate::TableCollection)
/// * `is_sample`: A `0/1` array indicating if each node in `tables`
///    is a sample (`1`) or is not a sample (`0`).
/// * `convert_time`: A callback to convert time, *e.g.* from forwards
///                   to backwards. For example, see
///                   [``simple_time_reverser``](crate::tskit_tools::simple_time_reverser).
/// * `build_indexes`: If `true`, build the edge table indexes for the return value.
///
/// # Notes
///
/// If the input ``tables`` are not sorted, pass ``false`` for
/// `build_indexes`.
///
/// This function will not be part of the long-term API.
/// Rather, it is the minimum currently needed to get stuff done.
///
/// # Returns
///
/// A [``tskit::TableCollection``].
///
/// # Example
///
/// ```
/// use tskit::TableAccess;
/// let mut tables = forrustts::TableCollection::new(100).unwrap();
/// tables.add_node(0., 0).unwrap(); // Add a parent node at time 0
/// tables.add_node(1., 0).unwrap(); // Add a child node at time 1
/// tables.add_edge(0, 100, 0, 1).unwrap(); // Add an edge
/// let is_sample = vec![0, 1]; // Mark the child node as a sample.
/// let tsk_tables = forrustts::tskit_tools::convert_to_tskit_minimal(
///     &tables,
///     &is_sample,
///     forrustts::tskit_tools::simple_time_reverser(1.),
///     true,
/// );
/// assert_eq!(tsk_tables.nodes().num_rows(), 2);
/// assert_eq!(tsk_tables.edges().num_rows(), 1);
/// assert_eq!(tsk_tables.populations().num_rows(), 1);
/// ```
pub fn convert_to_tskit_minimal(
    tables: &TableCollection,
    is_sample: &[i32],
    convert_time: impl Fn(Time) -> f64,
    build_indexes: bool,
) -> tskit::TableCollection {
    let mut tsk_tables = tskit::TableCollection::new(tables.genome_length().0 as f64).unwrap();

    for e in tables.edges() {
        tsk_tables
            .add_edge(e.left.0 as f64, e.right.0 as f64, e.parent.0, e.child.0)
            .unwrap();
    }

    let mut max_pop: tsk_id_t = -1;
    for (i, n) in tables.enumerate_nodes() {
        let flags: tsk_flags_t = if is_sample[i] > 0 {
            TSK_NODE_IS_SAMPLE
        } else {
            0
        };
        tsk_tables
            .add_node(flags, convert_time(n.time), n.deme.0, TSK_NULL)
            .unwrap();
        max_pop = std::cmp::max(n.deme.0, max_pop);
    }

    for _ in 0..(max_pop + 1) {
        tsk_tables.add_population().unwrap();
    }

    if build_indexes {
        tsk_tables.build_index().unwrap();
    }

    tsk_tables
}

fn swap_with_empty<T>(v: &mut Vec<T>) {
    let mut temp = Vec::<T>::new();
    std::mem::swap(v, &mut temp);
}

/// Convert a [``TableCollection``](crate::TableCollection)
/// to ``tskit`` format.
///
/// The input tables are mutable and have all of their
/// data "cleaned out".
///
/// # Parameters
///
/// * `is_sample`: A `0/1` array indicating if each node in `tables`
///    is a sample (`1`) or is not a sample (`0`).
/// * `convert_time`: A callback to convert time, *e.g.* from forwards
///                   to backwards. For example, see
///                   [``simple_time_reverser``](crate::tskit_tools::simple_time_reverser).
/// * `build_indexes`: If `true`, build the edge table indexes for the return value.
/// * `tables`: A mutable [``TableCollection``](crate::TableCollection)
///
/// # Notes
///
/// If the input ``tables`` are not sorted, pass ``false`` for
/// `build_indexes`.
///
/// This function will not be part of the long-term API.
/// Rather, it is the minimum currently needed to get stuff done.
///
/// # Returns
///
/// A [``tskit::TableCollection``].
///
/// # Example
///
/// ```
/// use tskit::TableAccess;
/// let mut tables = forrustts::TableCollection::new(100).unwrap();
/// tables.add_node(0., 0).unwrap(); // Add a parent node at time 0
/// tables.add_node(1., 0).unwrap(); // Add a child node at time 1
/// tables.add_edge(0, 100, 0, 1).unwrap(); // Add an edge
/// let is_sample = vec![0, 1]; // Mark the child node as a sample.
/// let tsk_tables = forrustts::tskit_tools::convert_to_tskit_and_drain_minimal(
///     &is_sample,
///     forrustts::tskit_tools::simple_time_reverser(1.),
///     true,
///     &mut tables,
/// );
/// assert_eq!(tsk_tables.nodes().num_rows(), 2);
/// assert_eq!(tsk_tables.edges().num_rows(), 1);
/// assert_eq!(tsk_tables.populations().num_rows(), 1);
/// // The input tables have no data:
/// assert_eq!(tables.num_nodes(), 0);
/// assert_eq!(tables.num_edges(), 0);
/// ```
pub fn convert_to_tskit_and_drain_minimal(
    is_sample: &[i32],
    convert_time: impl Fn(Time) -> f64,
    build_indexes: bool,
    tables: &mut TableCollection,
) -> tskit::TableCollection {
    let mut tsk_tables = tskit::TableCollection::new(tables.genome_length().0 as f64).unwrap();

    let mut max_pop: tsk_id_t = -1;
    for (i, n) in tables.enumerate_nodes() {
        let flags: tsk_flags_t = if is_sample[i] > 0 {
            TSK_NODE_IS_SAMPLE
        } else {
            0
        };
        tsk_tables
            .add_node(flags, convert_time(n.time), n.deme.0, TSK_NULL)
            .unwrap();
        max_pop = std::cmp::max(n.deme.0, max_pop);
    }
    swap_with_empty(&mut tables.nodes_);

    // Edges take the most memory,
    // so we clear out the other tables first.
    for e in tables.edges() {
        tsk_tables
            .add_edge(e.left.0 as f64, e.right.0 as f64, e.parent.0, e.child.0)
            .unwrap();
    }
    swap_with_empty(&mut tables.edges_);

    for _ in 0..(max_pop + 1) {
        tsk_tables.add_population().unwrap();
    }

    if build_indexes {
        tsk_tables.build_index().unwrap();
    }

    assert_eq!(tables.nodes_.capacity(), 0);
    assert_eq!(tables.edges_.capacity(), 0);

    tsk_tables
}

#[cfg(test)]
mod tests {

    use super::*;
    use tskit::TableAccess;

    #[test]
    fn test_convert_to_tskit() {
        let mut tables = TableCollection::new(100).unwrap();
        tables.add_node(0., 0).unwrap(); // Add a parent node at time 0
        tables.add_node(1., 0).unwrap(); // Add a child node at time 1
        tables.add_edge(0, 100, 0, 1).unwrap(); // Add an edge
        let is_sample = vec![0, 1]; // Mark the child node as a sample.
        let tsk_tables =
            convert_to_tskit_minimal(&tables, &is_sample, simple_time_reverser(1.), true);
        assert_eq!(tsk_tables.nodes().num_rows(), 2);
        assert_eq!(tsk_tables.edges().num_rows(), 1);
        assert_eq!(tsk_tables.populations().num_rows(), 1);
        assert_eq!(tables.num_edges(), 1);
        assert_eq!(tables.num_nodes(), 2);
    }

    #[test]
    fn test_convert_to_tskit_and_drain() {
        let mut tables = TableCollection::new(100).unwrap();
        tables.add_node(0., 0).unwrap(); // Add a parent node at time 0
        tables.add_node(1., 0).unwrap(); // Add a child node at time 1
        tables.add_edge(0, 100, 0, 1).unwrap(); // Add an edge
        let is_sample = vec![0, 1]; // Mark the child node as a sample.
        let tsk_tables = convert_to_tskit_and_drain_minimal(
            &is_sample,
            simple_time_reverser(1.),
            true,
            &mut tables,
        );
        assert_eq!(tsk_tables.nodes().num_rows(), 2);
        assert_eq!(tsk_tables.edges().num_rows(), 1);
        assert_eq!(tsk_tables.populations().num_rows(), 1);
        assert_eq!(tables.edges_.capacity(), 0);
        assert_eq!(tables.nodes_.capacity(), 0);
    }
}
