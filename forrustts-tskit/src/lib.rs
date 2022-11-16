#![warn(missing_docs)]

//! Data interchange with ``tskit`` format using
//! [``tskit``](https://crates.io/crates/tskit).
//!
//! This crate provides low-level functions to convert
//! [`forrustts_tables_trees::TableCollection`] to
//! [`TskTableCollection`].
//!
//! For simple cases, [`export_tables`] will suffice.
//! But, when metadata export is needed, use the
//! per-table functions.
//!
//! In general, the best way to validate your output
//! will be to create a tree sequence from the
//! [`TskTableCollection`].  Doing so triggers
//! validation code.  See the examples [here](export_tables).
//!
//! # Technical notes
//!
//! When needed, metadata must be stored in containers
//! that are accessible via `usize` via the trait
//! [`GetWithUsize`].  That trait is implemented
//! for `Vec` and `HashMap`.
//!
//! ## Limitations
//!
//! * We do not have a good model for going back and forth
//!   between nodes and individuals.
//! * Possible logic errors are probably untested.
//! * It would be preferable to do this with procedural
//!   macros if at all possible.

use bitflags::bitflags;
use forrustts_core::newtypes::Time;
use forrustts_core::traits::TableType;
use forrustts_tables_trees::TableCollection;
use thiserror::Error;

fn swap_with_empty<T>(v: &mut Vec<T>) {
    let mut temp = Vec::<T>::new();
    std::mem::swap(v, &mut temp);
}

/// Type alias for [`tskit::TableCollection`].
/// This exists to make the documentation less ambiguous.
pub type TskTableCollection = tskit::TableCollection;

/// In order to support a wider range of
/// metadata containers, we provide a trait
/// that is an abstraction over the `get` functions
/// found in `std` containers.
///
/// We require `usize` because we will need to
/// match up the `i-th` value in a table
/// with its metadata, and table indexes/enumerations
/// use `usize`.
///
/// This trait allows one to use either `Vec` or
/// `HashMap` to store metadata.
pub trait GetWithUsize {
    /// The type to return
    type Output: Sized;
    /// Get an optional reference to `Self::Output`.
    /// The implementation should have the behavior
    /// of the `get` functions of `Vec` and/or `HashMap`.
    fn get(&self, index: usize) -> Option<&Self::Output>;
}

impl<T> GetWithUsize for &[T] {
    type Output = T;
    fn get(&self, index: usize) -> Option<&Self::Output> {
        <[T]>::get(self, index)
    }
}

impl<T> GetWithUsize for Vec<T> {
    type Output = T;
    fn get(&self, index: usize) -> Option<&Self::Output> {
        self.as_slice().get(index)
    }
}

impl<V, S> GetWithUsize for std::collections::HashMap<usize, V, S>
where
    S: std::hash::BuildHasher,
{
    type Output = V;
    fn get(&self, index: usize) -> Option<&Self::Output> {
        std::collections::HashMap::<usize, V, S>::get(self, &index)
    }
}

bitflags! {
    /// Flags affecting the behavior of [`export_tables`].
    ///
    /// The default/empty state is to "do nothing".
    #[derive(Default)]
    pub struct TableCollectionExportFlags: u32 {
        /// Build the edge table indexes
        const BUILD_INDEXES = 1 << 0;
    }
}

/// Error type returned by [``export_tables``]
#[derive(Error, Debug)]
pub enum TableCollectionExportError {
    /// Error returned from `tskit`
    #[error("{value:?}")]
    TskitError {
        /// The specific error returned from
        /// `tskit`.
        #[from]
        value: tskit::TskitError,
    },
    /// Returned when a metadata key cannot fetch
    /// a metadata object.
    #[error("Invalid metadata key")]
    InvalidMetadataKey,
    /// Returned when a table is unexpectedly empty.
    #[error("Empty table")]
    EmptyTable,
}

/// Return a closure to help reverse time.
///
/// For all input values, ``t`` the closure will
/// return ``-1.0*(t - x) as f64``.
pub fn simple_time_reverser<T: Into<Time>>(x: T) -> impl Fn(Time) -> f64 {
    let ti = x.into();
    move |t: Time| -1. * (t.into_raw() - ti.into_raw()) as f64
}

/// Convert a [``TableCollection``](crate::TableCollection)
/// to ``tskit`` format.
///
/// # Parameters
///
/// * `tables`: A [``TableCollection``](crate::TableCollection)
/// * `convert_time`: A callback to convert time, *e.g.* from forwards
///                   to backwards. For example, see
///                   [``simple_time_reverser``](simple_time_reverser).
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
/// A [``TskTableCollection``].
///
/// # Example
///
/// ```
/// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// // Add 2 sample nodes
/// tables.add_node_with_flags(0., 0,
///                            forrustts_tables_trees::NodeFlags::IS_SAMPLE.bits()).unwrap();
/// tables.add_node_with_flags(1., 0,
///                            forrustts_tables_trees::NodeFlags::IS_SAMPLE.bits()).unwrap();
/// tables.add_edge(0, 100, 0, 1).unwrap(); // Add an edge
/// let tsk_tables = forrustts_tskit::export_tables(
///     tables,
///     &forrustts_tskit::simple_time_reverser(1),
///     forrustts_tskit::TableCollectionExportFlags::BUILD_INDEXES,
/// ).unwrap();
/// assert_eq!(tsk_tables.nodes().num_rows(), 2);
/// assert_eq!(tsk_tables.edges().num_rows(), 1);
/// assert_eq!(tsk_tables.populations().num_rows(), 1);
///
/// // Creating a tree sequence is the best means
/// // of validating.
/// let ts = tsk_tables.tree_sequence(tskit::TreeSequenceFlags::default()).unwrap();
/// ```
pub fn export_tables<F: Into<Option<TableCollectionExportFlags>>>(
    tables: TableCollection,
    convert_time: impl Fn(Time) -> f64,
    export_flags: F,
) -> Result<TskTableCollection, TableCollectionExportError> {
    let flags = match export_flags.into() {
        Some(x) => x,
        None => TableCollectionExportFlags::default(),
    };

    let mut tables_copy = tables;
    let mut tsk_tables = TskTableCollection::new(tables_copy.genome_length().into_raw() as f64)?;

    export_edges(tables_copy.edges(), &mut tsk_tables)?;
    let mut edges = tables_copy.dump_edge_table();
    swap_with_empty(&mut edges);
    debug_assert_eq!(tables_copy.edges().len(), 0);

    build_population_table(tables_copy.nodes(), &mut tsk_tables)?;
    export_nodes(tables_copy.nodes(), &convert_time, &mut tsk_tables)?;
    let mut nodes = tables_copy.dump_node_table();
    swap_with_empty(&mut nodes);
    debug_assert_eq!(tables_copy.nodes().len(), 0);

    export_sites(tables_copy.sites(), &mut tsk_tables)?;
    export_mutations(tables_copy.mutations(), &convert_time, &mut tsk_tables)?;

    if flags.contains(TableCollectionExportFlags::BUILD_INDEXES) {
        tsk_tables.build_index()?;
    }

    Ok(tsk_tables)
}

/// Export an edge table.
///
/// # Parameters
///
/// * `edges`: a slice of [`forrustts_tables_trees::Edge`]
/// * `tsk_tables`: The output tables.
///
/// # Example
///
/// ```
/// use forrustts_tables_trees::TableCollection;
/// use forrustts_core::traits::TableType;
/// use forrustts_tskit::TskTableCollection;
///
/// let mut tables = TableCollection::new(100).unwrap();
/// tables.add_edge(25, 50, 0, 1).unwrap();
/// let mut tsk_tables = TskTableCollection::new(tables.genome_length().into_raw() as f64).unwrap();
/// forrustts_tskit::export_edges(tables.edges(),
///                               &mut tsk_tables).unwrap();
/// assert_eq!(tsk_tables.edges().parent(0).unwrap(), 0);
/// assert_eq!(tsk_tables.edges().child(0).unwrap(), 1);
/// assert_eq!(tsk_tables.edges().left(0).unwrap(), 25.);
/// assert_eq!(tsk_tables.edges().right(0).unwrap(), 50.);
/// ```
///
/// # Errors
///
/// [`tskit::TskitError`] if adding rows returns an error.
pub fn export_edges(
    edges: &[forrustts_tables_trees::Edge],
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError> {
    for edge in edges {
        tsk_tables.add_edge(
            edge.left.into_raw() as f64,
            edge.right.into_raw() as f64,
            edge.parent.into_raw(),
            edge.child.into_raw(),
        )?;
    }
    Ok(())
}

/// Build a population table.
///
/// The maximal `deme` value is obtained from the
/// input node table.
/// The output tables are populated with enough rows
/// such that all deme IDs <= the maximal value
/// are valid.
///
/// # Parameters
///
/// * `nodes`: a slice of [`forrustts_tables_trees::Node`]
/// * `tsk_tables`: The output tables.
///
/// # Example
///
/// ```
/// use forrustts_tables_trees::TableCollection;
/// use forrustts_core::traits::TableType;
/// use forrustts_tskit::TskTableCollection;
///
/// let mut tables = TableCollection::new(100).unwrap();
/// tables.add_node(1, 2).unwrap();
/// let mut tsk_tables = TskTableCollection::new(tables.genome_length().into_raw() as f64).unwrap();
/// forrustts_tskit::build_population_table(tables.nodes(),
///                                         &mut tsk_tables).unwrap();
/// // Max deme ID is 2, so there should be 3 rows.
/// assert_eq!(tsk_tables.populations().num_rows(), 3);
/// ```
///
/// # Errors
///
/// [`TableCollectionExportError::EmptyTable`] if `nodes` is empty.
/// [`tskit::TskitError`] if adding rows returns an error.
///
/// ## Example error
///
/// ```should_panic
/// let tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// let mut tsk_tables = tskit::TableCollection::new(100.).unwrap();
/// forrustts_tskit::build_population_table(tables.nodes(), &mut tsk_tables).unwrap();
/// ```
pub fn build_population_table(
    nodes: &[forrustts_tables_trees::Node],
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError> {
    if nodes.is_empty() {
        return Err(TableCollectionExportError::EmptyTable);
    }
    let mut max_deme = forrustts_core::newtypes::DemeId::NULL;
    for node in nodes {
        max_deme = std::cmp::max(node.deme, max_deme);
    }
    for _ in 0..(max_deme.into_raw() + 1) {
        tsk_tables.add_population()?;
    }

    Ok(())
}

/// Builds a population table with metadata.
///
/// # Parameters
///
/// * `nodes`: A slice of [`forrustts_tables_trees::NodeId`]
/// * `metadata`: metadata container.
/// * `tsk_tables`: the output tables.
///
/// # Examples
///
/// ```
/// #[derive(serde::Serialize, serde::Deserialize, tskit::metadata::PopulationMetadata)]
/// #[serializer("serde_json")]
/// struct PopulationName(String);
///
/// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// tables.add_node(0, 0).unwrap();
/// let md = vec![PopulationName("YRB".to_string())];
///
/// let mut tsk_tables = tskit::TableCollection::new(100.).unwrap();
/// forrustts_tskit::build_population_table_with_metadata(tables.nodes(), &md, &mut tsk_tables).unwrap();
/// let decoded = tsk_tables.populations().metadata::<PopulationName>(0.into()).unwrap().unwrap();
/// assert_eq!(&decoded.0, "YRB");
/// ```
///
/// # Errors
///
/// [`TableCollectionExportError`] if `nodes` is empty, if
/// a deme does not have associated metadata,
/// or if there is an error from `tskit`.
///
pub fn build_population_table_with_metadata<'metadata, C, M>(
    nodes: &[forrustts_tables_trees::Node],
    metadata: &'metadata C,
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError>
where
    M: tskit::metadata::PopulationMetadata + 'metadata + Sized,
    C: GetWithUsize<Output = M>,
{
    if nodes.is_empty() {
        return Err(TableCollectionExportError::EmptyTable);
    }
    let mut max_deme = forrustts_core::newtypes::DemeId::NULL;
    for node in nodes {
        max_deme = std::cmp::max(node.deme, max_deme);
    }
    for p in 0..(max_deme.into_raw() + 1) {
        match metadata.get(p as usize) {
            Some(md) => tsk_tables.add_population_with_metadata(md)?,
            None => return Err(TableCollectionExportError::InvalidMetadataKey),
        };
    }

    Ok(())
}

// TODO: handle "individuals"

/// Export a node table.
///
/// # Parameters
///
/// * `nodes`: a slice of [`forrustts_tables_trees::Node`]
/// * `convert_time`: A function to convert input time (forwards) to `tskit`-time (backwards).
/// * `tsk_tables`: The output tables.
///
/// # Example
///
/// ```
/// use forrustts_tables_trees::TableCollection;
/// use forrustts_core::traits::TableType;
/// use forrustts_tskit::TskTableCollection;
///
/// let mut tables = TableCollection::new(100).unwrap();
/// tables.add_node(1, 2).unwrap();
/// let mut tsk_tables = TskTableCollection::new(tables.genome_length().into_raw() as f64).unwrap();
/// forrustts_tskit::export_nodes(tables.nodes(),
///                               &forrustts_tskit::simple_time_reverser(1),
///                               &mut tsk_tables).unwrap();
/// // Time 1 is our maximal time in the forward direction,
/// // so it is converted to 0:
/// assert_eq!(tsk_tables.nodes().time(0).unwrap(), 0.0);
/// assert_eq!(tsk_tables.nodes().population(0).unwrap(), 2);
/// ```
///
/// # Errors
///
/// [`tskit::TskitError`] if adding rows returns an error.
pub fn export_nodes(
    nodes: &[forrustts_tables_trees::Node],
    convert_time: &impl Fn(Time) -> f64,
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError> {
    for node in nodes {
        tsk_tables.add_node(
            node.flags,
            convert_time(node.time),
            node.deme.into_raw(),
            tskit::IndividualId::NULL,
        )?;
    }
    Ok(())
}

// TODO: handle "individuals"
/// Export node table with metadata.
///
/// # Parameters
///
/// * `nodes`: a slice of [`forrustts_tables_trees::Node`]
/// * `convert_time`: A function to convert input time (forwards) to `tskit`-time (backwards).
/// * `metadata`: the node metadata
/// * `tsk_tables`: The output tables.
///
/// # Example
///
/// ```
/// use forrustts_tables_trees::TableCollection;
/// use forrustts_core::traits::TableType;
/// use forrustts_tskit::TskTableCollection;
/// use serde::{Serialize, Deserialize};
/// use tskit::metadata::NodeMetadata;
///
/// #[derive(Serialize, Deserialize, NodeMetadata)]
/// #[serializer("serde_json")]
/// struct Metadata(i32);
///
/// let mut tables = TableCollection::new(100).unwrap();
/// tables.add_node(1, 2).unwrap();
/// let mut tsk_tables = TskTableCollection::new(tables.genome_length().into_raw() as f64).unwrap();
/// let md = vec![Metadata(42)];
/// forrustts_tskit::export_nodes_with_metadata(tables.nodes(),
///                               &forrustts_tskit::simple_time_reverser(1),
///                               &md,
///                               &mut tsk_tables).unwrap();
/// // Time 1 is our maximal time in the forward direction,
/// // so it is converted to 0:
/// assert_eq!(tsk_tables.nodes().time(0).unwrap(), 0.0);
/// assert_eq!(tsk_tables.nodes().population(0).unwrap(), 2);
/// ```
///
/// # Errors
///
/// [`tskit::TskitError`] if adding rows returns an error.
pub fn export_nodes_with_metadata<'metadata, C, M>(
    nodes: &[forrustts_tables_trees::Node],
    convert_time: &impl Fn(Time) -> f64,
    metadata: &'metadata C,
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError>
where
    M: tskit::metadata::NodeMetadata + 'metadata + Sized,
    C: GetWithUsize<Output = M>,
{
    for (i, node) in nodes.iter().enumerate() {
        match metadata.get(i) {
            Some(md) => tsk_tables.add_node_with_metadata(
                node.flags,
                convert_time(node.time),
                node.deme.into_raw(),
                tskit::IndividualId::NULL,
                md,
            )?,
            None => return Err(TableCollectionExportError::InvalidMetadataKey),
        };
    }
    Ok(())
}

/// Export a mutation table
///
/// # Examples
///
/// ```
/// use forrustts_tskit::simple_time_reverser;
///
/// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// tables.add_mutation(0, None, 0, 0, None, true).unwrap();
/// let mut tsk_tables = tskit::TableCollection::new(100.).unwrap();
/// forrustts_tskit::export_mutations(tables.mutations(), &simple_time_reverser(0), &mut tsk_tables).unwrap();
/// assert_eq!(tsk_tables.mutations().num_rows(), 1);
/// ```
pub fn export_mutations(
    mutations: &[forrustts_tables_trees::MutationRecord],
    convert_time: &impl Fn(Time) -> f64,
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError> {
    for mutation in mutations {
        tsk_tables.add_mutation(
            mutation.site.into_raw(),
            mutation.node.into_raw(),
            tskit::MutationId::NULL,
            convert_time(mutation.time),
            match &mutation.derived_state {
                Some(x) => Some(x),
                None => None,
            },
        )?;
    }
    Ok(())
}

/// Export a mutation table with metadata
///
/// # Examples
///
/// ```
/// use serde::{Serialize, Deserialize};
/// use tskit::metadata::MutationMetadata;
/// use forrustts_tskit::simple_time_reverser;
///
/// #[derive(Serialize, Deserialize, MutationMetadata)]
/// #[serializer("serde_json")]
/// struct Metadata(i32);
///
///
/// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// tables.add_mutation(0, None, 0, 0, None, true).unwrap();
/// let mut tsk_tables = tskit::TableCollection::new(100.).unwrap();
/// let md = vec![Metadata(-1)];
/// forrustts_tskit::export_mutations_with_metadata(tables.mutations(), &simple_time_reverser(0), &md, &mut tsk_tables).unwrap();
/// assert_eq!(tsk_tables.mutations().num_rows(), 1);
/// ```
pub fn export_mutations_with_metadata<'metadata, C, M>(
    mutations: &[forrustts_tables_trees::MutationRecord],
    convert_time: &impl Fn(Time) -> f64,
    metadata: &'metadata C,
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError>
where
    M: tskit::metadata::MutationMetadata + 'metadata + Sized,
    C: GetWithUsize<Output = M>,
{
    for mutation in mutations {
        match mutation.key {
            Some(key) => match metadata.get(key) {
                Some(md) => {
                    tsk_tables.add_mutation_with_metadata(
                        mutation.site.into_raw(),
                        mutation.node.into_raw(),
                        tskit::MutationId::NULL,
                        convert_time(mutation.time),
                        match &mutation.derived_state {
                            Some(x) => Some(x),
                            None => None,
                        },
                        md,
                    )?;
                }
                None => {
                    return Err(TableCollectionExportError::InvalidMetadataKey);
                }
            },
            None => {
                tsk_tables.add_mutation(
                    mutation.site.into_raw(),
                    mutation.node.into_raw(),
                    tskit::MutationId::NULL,
                    mutation.time.into_raw() as f64,
                    match &mutation.derived_state {
                        Some(x) => Some(x),
                        None => None,
                    },
                )?;
            }
        }
    }
    Ok(())
}

/// Export a site table
///
/// # Examples
///
/// ```
/// use forrustts_tskit::simple_time_reverser;
///
/// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// tables.add_site(50, None).unwrap();
/// let mut tsk_tables = tskit::TableCollection::new(100.).unwrap();
/// forrustts_tskit::export_sites(tables.sites(), &mut tsk_tables).unwrap();
/// assert_eq!(tsk_tables.sites().num_rows(), 1);
/// ```
pub fn export_sites(
    sites: &[forrustts_tables_trees::Site],
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError> {
    for site in sites {
        tsk_tables.add_site(
            site.position.into_raw() as f64,
            match &site.ancestral_state {
                Some(x) => Some(x),
                None => None,
            },
        )?;
    }
    Ok(())
}

/// Export a site table
///
/// # Examples
///
/// ```
/// use serde::{Serialize, Deserialize};
/// use tskit::metadata::SiteMetadata;
/// use forrustts_tskit::simple_time_reverser;
///
/// #[derive(Serialize, Deserialize, SiteMetadata)]
/// #[serializer("serde_json")]
/// struct Metadata(i32);
///
/// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// tables.add_site(50, None).unwrap();
/// let mut tsk_tables = tskit::TableCollection::new(100.).unwrap();
/// let md = vec![Metadata(1234)];
/// forrustts_tskit::export_sites_with_metadata(tables.sites(), &md, &mut tsk_tables).unwrap();
/// assert_eq!(tsk_tables.sites().num_rows(), 1);
/// ```
pub fn export_sites_with_metadata<'metadata, C, M>(
    sites: &[forrustts_tables_trees::Site],
    metadata: &'metadata C,
    tsk_tables: &mut TskTableCollection,
) -> Result<(), TableCollectionExportError>
where
    M: tskit::metadata::SiteMetadata + 'metadata + Sized,
    C: GetWithUsize<Output = M>,
{
    for (i, site) in sites.iter().enumerate() {
        match metadata.get(i) {
            Some(x) => tsk_tables.add_site_with_metadata(
                site.position.into_raw() as f64,
                match &site.ancestral_state {
                    Some(x) => Some(x),
                    None => None,
                },
                x,
            )?,
            None => return Err(TableCollectionExportError::InvalidMetadataKey),
        };
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_node_flags_against_tskit() {
        // We cannot export our flags directly
        // unless we agree with tskit on some definitions
        assert_eq!(
            forrustts_tables_trees::NodeFlags::IS_SAMPLE.bits(),
            tskit::TSK_NODE_IS_SAMPLE
        );
    }
}
