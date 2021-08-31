use crate::newtypes::{
    DemeId, EdgeId, MutationId, NodeId, Position, SiteId, TablesIdInteger, Time,
};
use bitflags::bitflags;
use std::cmp::Ordering;
use thiserror::Error;

/// Error type related to [``TableCollection``]
#[derive(Error, Debug, PartialEq)]
pub enum TablesError {
    /// Returned by [``TableCollection::new``].
    #[error("Invalid genome length")]
    InvalidGenomeLength,
    /// Returned when invalid node `ID`s are encountered.
    #[error("Invalid node: {found:?}")]
    InvalidNodeValue {
        /// The invalid `ID`
        found: NodeId,
    },
    /// Returned when invalid positions are encountered.
    #[error("Invalid value for position: {found:?}")]
    InvalidPosition {
        /// The invalid position
        found: Position,
    },
    /// Returned when table validation detects duplicate posititions
    /// in a site table.
    #[error("Duplicated site positions found")]
    DuplicatedSitePosition,
    /// Returned then site tables are not properly sorted
    #[error("Site positions are unsorted")]
    UnsortedSitePosition,
    #[error("Site ID out of bounds")]
    /// Returned when a [``MutationRecord``]'s [`SiteId`] is out of bounds.
    SiteOutofBounds,
    /// Returned when mutations tables are not sorted by site position.
    #[error("Mutations not sorted by increasing position")]
    UnsortedMutationPositions,
    /// Returned when mutations at the same site are not sorted
    /// correctly by time.
    #[error("Mutations within same site are not sorted by time")]
    UnsortedMutationsWithinSite,
    /// Retured when a [``MutationRecord``]'s time field is not finite
    #[error("Invalid Mutation time.")]
    InvalidMutationTime,
    /// Retured when a [``Node``]'s time field is not finite,
    /// including the node field of [``MutationRecord``].
    #[error("Invalid Node time.")]
    InvalidNodeTime,
    /// Returned when an [``Edge``]'s left/right
    /// values are invalid.
    #[error("Invalid position range: {found:?}")]
    InvalidLeftRight {
        /// The invalid `(left, right)`.
        found: (Position, Position),
    },
    /// Returned when invalid times are encountered.
    #[error("Invalid value for time: {found:?}")]
    InvalidTime {
        /// The invalid time
        found: Time,
    },
    #[error("Invalid value for deme: {found:?}")]
    /// Returned with a deme's `ID` is invalid.
    InvalidDeme {
        /// The invalide deme `ID`
        found: DemeId,
    },
    #[error("Parent is NULL_ID")]
    /// Can be returned by [``validate_edge_table``]
    NullParent,
    #[error("Child is NULL_ID")]
    /// Can be returned by [``validate_edge_table``]
    NullChild,
    #[error("Node is out of bounds")]
    /// Can be returned by [``validate_edge_table``]
    NodeOutOfBounds,
    #[error("Node time order violation")]
    /// Can be returned by [``validate_edge_table``]
    NodeTimesUnordered,
    #[error("Parents not sorted by time")]
    /// Can be returned by [``validate_edge_table``]
    ParentTimesUnsorted,
    #[error("Parents not contiguous")]
    /// Can be returned by [``validate_edge_table``]
    ParentsNotContiguous,
    #[error("Edges not sorted by child")]
    /// Can be returned by [``validate_edge_table``]
    EdgesNotSortedByChild,
    #[error("Edges not sorted by left")]
    /// Can be returned by [``validate_edge_table``]
    EdgesNotSortedByLeft,
    #[error("Duplicate edges")]
    /// Can be returned by [``validate_edge_table``]
    DuplicateEdges,
    /// Can be returned by [`crate::TreeSequence::new`]
    /// and [`crate::TreeSequence::new_with_samples`]
    #[error("Tables not indexed")]
    TablesNotIndexed,
}

/// Result type for operations on tables
pub type TablesResult<T> = std::result::Result<T, TablesError>;

/// A Node of a tree sequence
#[derive(Copy, Clone)]
pub struct Node {
    /// Birth time
    pub time: Time,
    /// Population (deme) of node
    pub deme: DemeId,
    /// Bit flags
    pub flags: u32,
}

/// An Edge is a transmission event
///
/// An edge is a record of transmission of
/// a half-open chunk of genome `[left, right)`
/// from `parent` to `child`.
#[derive(Copy, Clone)]
pub struct Edge {
    /// Left end
    pub left: Position,
    /// Right end
    pub right: Position,
    /// Index of parent in a [NodeTable](type.NodeTable.html)
    pub parent: NodeId,
    /// Index of child in a [NodeTable](type.NodeTable.html)
    pub child: NodeId,
}

// TODO: It would be nice to use generics here
// to allow arbitrary types for ancestral_state
// and derived_state.

/// A Site is the location and
/// ancestral state of a tables::MutationRecord
#[derive(Clone)]
pub struct Site {
    /// Position of the mutation
    pub position: Position,
    /// The ancestral state.
    /// [``None``] implies client code
    /// will apply a default.
    pub ancestral_state: Option<Vec<u8>>,
}

/// A MutationRecord is the minimal information
/// needed about a mutation to track it
/// on a tree sequence.
#[derive(Clone)]
pub struct MutationRecord {
    /// The node where the mutation maps
    pub node: NodeId,
    /// Reference to mutation metadata.
    pub key: Option<usize>,
    /// The index of the corresponding [``Site``].
    pub site: SiteId,
    /// The origin time of the mutation
    pub time: Time,
    /// The derived state.
    /// [``None``] implies client code
    /// will apply a default.
    pub derived_state: Option<Vec<u8>>,
    /// [``true``] if mutation affects fitness,
    /// [``false``] otherwise.
    pub neutral: bool,
}

/// A node table
pub type NodeTable = Vec<Node>;
/// An edge table
pub type EdgeTable = Vec<Edge>;
/// A site table
pub type SiteTable = Vec<Site>;
/// A Mutation table
pub type MutationTable = Vec<MutationRecord>;

fn position_non_negative(x: Position) -> TablesResult<()> {
    if x.0 < 0 {
        Err(TablesError::InvalidPosition { found: x })
    } else {
        Ok(())
    }
}

fn node_non_negative(x: NodeId) -> TablesResult<()> {
    if x < 0 {
        Err(TablesError::InvalidNodeValue { found: x })
    } else {
        Ok(())
    }
}

fn edge_table_add_row(
    edges: &mut EdgeTable,
    left: Position,
    right: Position,
    parent: NodeId,
    child: NodeId,
) -> TablesResult<EdgeId> {
    if right <= left {
        return Err(TablesError::InvalidLeftRight {
            found: (left, right),
        });
    }
    position_non_negative(left)?;
    position_non_negative(right)?;
    node_non_negative(parent)?;
    node_non_negative(child)?;

    edges.push(Edge {
        left,
        right,
        parent,
        child,
    });

    Ok(EdgeId::from(edges.len() - 1))
}

// NOTE: we allow negative times, in order to support "precapitation".
fn node_table_add_row(
    nodes: &mut NodeTable,
    time: Time,
    deme: DemeId,
    flags: u32,
) -> TablesResult<NodeId> {
    nodes.push(Node { time, deme, flags });

    Ok(NodeId::from(nodes.len() - 1))
}

fn site_table_add_row(
    sites: &mut SiteTable,
    position: Position,
    ancestral_state: Option<Vec<u8>>,
) -> TablesResult<SiteId> {
    position_non_negative(position)?;
    sites.push(Site {
        position,
        ancestral_state,
    });

    Ok(SiteId::from(sites.len() - 1))
}

fn mutation_table_add_row(
    mutations: &mut MutationTable,
    node: NodeId,
    key: Option<usize>,
    site: SiteId,
    time: Time,
    derived_state: Option<Vec<u8>>,
    neutral: bool,
) -> TablesResult<MutationId> {
    node_non_negative(node)?;
    mutations.push(MutationRecord {
        node,
        key,
        site,
        time,
        derived_state,
        neutral,
    });

    Ok(MutationId::from(mutations.len() - 1))
}

fn sort_edges(nodes: &[Node], edges: &mut [Edge]) {
    edges.sort_by(|a, b| {
        let aindex = a.parent.0 as usize;
        let bindex = b.parent.0 as usize;
        let ta = nodes[aindex].time;
        let tb = nodes[bindex].time;
        match ta.partial_cmp(&tb) {
            Some(std::cmp::Ordering::Equal) => {
                if a.parent == b.parent {
                    if a.child == b.child {
                        return a.left.cmp(&b.left);
                    }
                    a.child.cmp(&b.child)
                } else {
                    a.parent.cmp(&b.parent)
                }
            }
            Some(x) => x.reverse(),
            None => panic!("invalid parent times"),
        }
    });
}

fn record_site(sites: &[Site], mutation: &mut MutationRecord, new_site_table: &mut SiteTable) {
    let position = sites[mutation.site.0 as usize].position;
    if new_site_table.is_empty() || new_site_table[new_site_table.len() - 1].position != position {
        new_site_table.push(sites[mutation.site.0 as usize].clone());
    }

    mutation.site = SiteId((new_site_table.len() - 1) as TablesIdInteger);
}

fn sort_mutation_table(sites: &[Site], mutations: &mut [MutationRecord]) {
    mutations.sort_by(|a, b| {
        let pa = sites[a.site.0 as usize].position;
        let pb = sites[b.site.0 as usize].position;
        match pa.cmp(&pb) {
            std::cmp::Ordering::Equal => match a.time.partial_cmp(&b.time) {
                Some(x) => x,
                None => panic!("bad mutation times {} {}", a.time.0, b.time.0),
            },
            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
        }
    });
}

bitflags! {
    /// Set properties of a [`Node`].
    ///
    /// The first 16 bits are reserved for internal use.
    /// Client code is free to use the remaining bits
    /// as needed.
    #[derive(Default)]
    pub struct NodeFlags: u32 {
        /// Default
        const NONE = 0;
        /// The node is a sample node.
        const IS_SAMPLE = 1 << 0;
        /// The node is alive.
        /// Usually, this is set along with
        /// IS_SAMPLE in order to distinguish
        /// living individuals from, e.g.,
        /// ancient samples.
        const IS_ALIVE = 1 << 1;
    }
}

bitflags! {
    /// Modifies behavior of
    /// [``TableCollection::validate``]
    ///
    /// ```
    /// let f = forrustts_tables_trees::TableValidationFlags::default();
    /// assert_eq!(f.contains(forrustts_tables_trees::TableValidationFlags::VALIDATE_ALL), true);
    /// ```
    pub struct TableValidationFlags: u32 {
        /// Validate the edge table
        const VALIDATE_EDGES = 1<<0;
        /// Validate the site table
        const VALIDATE_SITES = 1<<1;
        /// Validate the mutation table
        const VALIDATE_MUTATIONS = 1<<2;
        /// Validate the node table
        const VALIDATE_NODES = 1<<3;
        /// Validate all tables.
        /// This is also the "default" value.
        const VALIDATE_ALL = Self::VALIDATE_EDGES.bits|Self::VALIDATE_MUTATIONS.bits|Self::VALIDATE_SITES.bits|Self::VALIDATE_NODES.bits;
    }
}

impl Default for TableValidationFlags {
    fn default() -> Self {
        TableValidationFlags::VALIDATE_ALL
    }
}

bitflags! {
    /// Modifies behavior of
    /// [``TableCollection::sort_tables``]
    ///
    /// ```
    /// let f = forrustts_tables_trees::TableSortingFlags::empty();
    /// assert_eq!(f.contains(forrustts_tables_trees::TableSortingFlags::SORT_ALL), true);
    /// ```
    #[derive(Default)]
    pub struct TableSortingFlags: u32 {
        /// Sort all tables.
        /// This is also the "default"/empty.
        const SORT_ALL = 0;
        /// Do not sort the edge table.
        const SKIP_EDGE_TABLE = 1 << 0;
    }
}

bitflags! {
    /// Modifies behavior of
    /// [``TableCollection::build_indexes``]
    #[derive(Default)]
    pub struct IndexTablesFlags: u32 {
        /// Default behavior
        const NONE = 0;
        /// Do not validate edge table
        const NO_VALIDATION = 1<<0;
    }
}

/// Perform a data integrity check on an [``EdgeTable``].
///
/// This checks, amongst other things, the sorting order
/// of the edges.
///
/// # Parameters
///
/// * `len`, the genome length of the tables.
///          Best obtained via [``TableCollection::genome_length``].
/// * `edges`, the [``EdgeTable``]
/// * `nodes`, the [``NodeTable``]
///
/// # Return
///
/// Returns ``Ok(true)`` if the tables pass all tests.
/// This return value allows this function to be used in
/// things like [``debug_assert``].
///
/// # Errors
///
/// Will return [``TablesError``] if the tables are not valid.
///
/// # Example
///
/// ```
/// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
/// // (do some stuff now...)
/// let rv = forrustts_tables_trees::validate_edge_table(tables.genome_length(),
///                                         &tables.edges(),
///                                         &tables.nodes()).unwrap();
/// assert_eq!(rv, true);
/// ```
pub fn validate_edge_table(len: Position, edges: &[Edge], nodes: &[Node]) -> TablesResult<bool> {
    if edges.is_empty() {
        return Ok(true);
    }
    let mut parent_seen = vec![0; nodes.len()];
    let mut last_parent: usize = edges[0].parent.0 as usize;
    let mut last_child: usize = edges[0].child.0 as usize;
    let mut last_left: Position = edges[0].left;

    for (i, edge) in edges.iter().enumerate() {
        if edge.parent == NodeId::NULL {
            return Err(TablesError::NullParent);
        }
        if edge.child == NodeId::NULL {
            return Err(TablesError::NullChild);
        }
        if edge.parent < 0 || edge.parent.0 as usize >= nodes.len() {
            return Err(TablesError::NodeOutOfBounds);
        }
        if edge.child < 0 || edge.child.0 as usize >= nodes.len() {
            return Err(TablesError::NodeOutOfBounds);
        }
        if edge.left.0 < 0 || edge.left > len {
            return Err(TablesError::InvalidPosition { found: edge.left });
        }
        if edge.right.0 < 0 || edge.right > len {
            return Err(TablesError::InvalidPosition { found: edge.right });
        }
        if edge.left >= edge.right {
            return Err(TablesError::InvalidLeftRight {
                found: (edge.left, edge.right),
            });
        }

        // child time must be > parent time b/c time goes forwards
        if nodes[edge.child.0 as usize].time <= nodes[edge.parent.0 as usize].time {
            return Err(TablesError::NodeTimesUnordered);
        }

        if parent_seen[edge.parent.0 as usize] == 1 {
            return Err(TablesError::ParentsNotContiguous);
        }

        if i > 0 {
            match nodes[edge.parent.0 as usize]
                .time
                .partial_cmp(&nodes[last_parent].time)
            {
                Some(std::cmp::Ordering::Greater) => {
                    return Err(TablesError::ParentTimesUnsorted);
                }
                Some(std::cmp::Ordering::Equal) => {
                    if edge.parent.0 as usize == last_parent {
                        if (edge.child.0 as usize) < last_child {
                            return Err(TablesError::EdgesNotSortedByChild);
                        }
                        if edge.child.0 as usize == last_child {
                            match edge.left.cmp(&last_left) {
                                Ordering::Greater => (),
                                Ordering::Equal => return Err(TablesError::DuplicateEdges),
                                Ordering::Less => return Err(TablesError::EdgesNotSortedByLeft),
                            }
                        }
                    } else {
                        parent_seen[last_parent] = 1;
                    }
                }
                Some(_) => (),
                None => panic!("invalid node times"),
            }
        }
        last_parent = edge.parent.0 as usize;
        last_child = edge.child.0 as usize;
        last_left = edge.left;
    }

    Ok(true)
}

pub fn validate_node_table(nodes: &[Node]) -> TablesResult<()> {
    for n in nodes {
        if !n.time.0.is_finite() {
            return Err(TablesError::InvalidNodeTime);
        }
    }
    Ok(())
}

pub fn validate_site_table(len: Position, sites: &[Site]) -> TablesResult<()> {
    for (i, site) in sites.iter().enumerate() {
        if site.position < 0 || site.position >= len {
            return Err(TablesError::InvalidPosition {
                found: site.position,
            });
        }
        if i > 0 {
            if sites[i - 1].position == site.position {
                return Err(TablesError::DuplicatedSitePosition);
            }
            if sites[i - 1].position > site.position {
                return Err(TablesError::UnsortedSitePosition);
            }
        }
    }
    Ok(())
}

pub fn validate_mutation_table(
    mutations: &[MutationRecord],
    sites: &[Site],
    nodes: &[Node],
) -> TablesResult<()> {
    let mut last_site: Option<SiteId> = None;
    let mut last_time = Time::MIN;
    for (i, mutation) in mutations.iter().enumerate() {
        if !mutation.time.0.is_finite() {
            return Err(TablesError::InvalidMutationTime);
        }
        if mutation.site < 0 || (mutation.site.0 as usize) >= sites.len() {
            return Err(TablesError::SiteOutofBounds);
        }
        if mutation.node < 0 || (mutation.node.0 as usize) >= nodes.len() {
            return Err(TablesError::NodeOutOfBounds);
        }
        if !nodes[mutation.node.0 as usize].time.0.is_finite() {
            return Err(TablesError::InvalidNodeTime);
        }
        if i > 0 {
            if mutations[i - 1].site > mutation.site {
                return Err(TablesError::UnsortedMutationPositions);
            }
            if last_site.is_some() && Some(mutation.site) == last_site && mutation.time < last_time
            {
                return Err(TablesError::UnsortedMutationsWithinSite);
            }
        }
        last_site = Some(mutation.site);
        last_time = mutation.time;
    }
    Ok(())
}

/// A collection of node, edge, site, and mutation tables.
#[derive(Clone)]
pub struct TableCollection {
    length_: Position, // Not visible outside of this module

    pub(crate) nodes_: NodeTable,
    pub(crate) edges_: EdgeTable,
    pub(crate) sites_: SiteTable,
    pub(crate) mutations_: MutationTable,
    pub(crate) edge_input_order: Vec<usize>,
    pub(crate) edge_output_order: Vec<usize>,
    pub(crate) is_indexed: bool,
}

impl TableCollection {
    /// Create a new instance.
    ///
    /// # Parameters
    ///
    /// * `genome_length`: the total genome length for the tables.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if `genome_length < 1`.
    pub fn new<P: Into<Position>>(genome_length: P) -> TablesResult<TableCollection> {
        let p = genome_length.into();
        if p.0 < 1 {
            return Err(TablesError::InvalidGenomeLength);
        }

        Ok(TableCollection {
            length_: p,
            nodes_: NodeTable::new(),
            edges_: EdgeTable::new(),
            sites_: SiteTable::new(),
            mutations_: MutationTable::new(),
            edge_input_order: vec![],
            edge_output_order: vec![],
            is_indexed: false,
        })
    }

    /// Add a [``Node``] to the [``NodeTable``]
    ///
    /// # Parameters
    ///
    /// * `time`, the birth time.
    /// * `deme`, the deme where the node is found.
    ///
    /// # Returns
    ///
    /// A [``NodeId``].
    ///
    /// The value will be [``NodeId::NULL``] in the case of
    /// overflow.
    ///
    /// # Side effects
    ///
    /// Adding a node invalidates current table indexes.
    /// Therefore, this function results in [`TableCollection::is_indexed`]
    /// returning false.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if `deme < 0`.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// let id = tables.add_node(1. , 0).unwrap();
    /// assert_eq!(id, 0);
    /// ```
    pub fn add_node<T: Into<Time>, D: Into<DemeId> + Copy>(
        &mut self,
        time: T,
        deme: D,
    ) -> TablesResult<NodeId> {
        self.add_node_with_flags(time, deme, NodeFlags::default().bits())
    }

    /// Add a [``Node``] to the [``NodeTable``] with flags set.
    ///
    /// # Parameters
    ///
    /// * `time`: the birth time.
    /// * `deme`: the deme where the node is found.
    /// * `flags`: node flags.  See [`NodeFlags`].
    ///
    /// # Returns
    ///
    /// A [``NodeId``].
    ///
    /// The value will be [``NodeId::NULL``] in the case of
    /// overflow.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// let id = tables.add_node_with_flags(1., 0,
    ///     (forrustts_tables_trees::NodeFlags::IS_ALIVE | forrustts_tables_trees::NodeFlags::IS_SAMPLE).bits()).unwrap();
    /// assert_eq!(id, 0);
    /// assert!(tables.node(0).flags & forrustts_tables_trees::NodeFlags::IS_ALIVE.bits() > 0);
    /// assert!(tables.node(0).flags & forrustts_tables_trees::NodeFlags::IS_SAMPLE.bits() > 0);
    ///
    /// ```
    ///
    /// ## Modifying node flags
    ///
    /// To modify flags, we recommend a dump/modify/set work flow,
    /// illustrated in the next example.
    ///
    /// The dump/set operations are constant time, moving the relevant vectors.
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// let id = tables.add_node_with_flags(1., 0,
    ///     (forrustts_tables_trees::NodeFlags::IS_ALIVE | forrustts_tables_trees::NodeFlags::IS_SAMPLE).bits()).unwrap();
    /// assert_eq!(id, 0);
    /// assert!(tables.node(0).flags & forrustts_tables_trees::NodeFlags::IS_ALIVE.bits() > 0);
    /// assert!(tables.node(0).flags & forrustts_tables_trees::NodeFlags::IS_SAMPLE.bits() > 0);
    ///
    /// let mut nodes = tables.dump_nodes();
    /// assert_eq!(tables.num_nodes(), 0);
    /// for n in &mut nodes {
    ///     n.flags = 0; // reset all the flags
    /// }
    /// tables.set_node_table(nodes);
    /// assert!(tables.node(0).flags & forrustts_tables_trees::NodeFlags::IS_ALIVE.bits() == 0);
    /// assert!(tables.node(0).flags & forrustts_tables_trees::NodeFlags::IS_SAMPLE.bits() == 0);
    /// ```
    pub fn add_node_with_flags<T: Into<Time>, D: Into<DemeId> + Copy>(
        &mut self,
        time: T,
        deme: D,
        flags: u32,
    ) -> TablesResult<NodeId> {
        self.is_indexed = false;
        node_table_add_row(&mut self.nodes_, time.into(), deme.into(), flags)
    }

    /// Add an [``Edge``] to the [``EdgeTable``].
    ///
    /// # Parameters
    ///
    /// * `left`, the left end of the edge
    /// * `right`, the right end of the edge
    /// * `parent`, the parent of the edge
    /// * `child`, the child of the edge
    ///
    /// # Returns
    ///
    /// An [``EdgeId``].
    ///
    /// The value will be [``EdgeId::NULL``] in the case of
    /// overflow.
    ///
    ///
    /// # Side effects
    ///
    /// Adding an edge invalidates current table indexes.
    /// Therefore, this function results in [`TableCollection::is_indexed`]
    /// returning false.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if any of the input
    /// are invalid.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// let id = tables.add_edge(0, 3, 5, 9).unwrap();
    /// assert_eq!(id, 0);
    /// ```
    pub fn add_edge<L: Into<Position>, R: Into<Position>, P: Into<NodeId>, C: Into<NodeId>>(
        &mut self,
        left: L,
        right: R,
        parent: P,
        child: C,
    ) -> TablesResult<EdgeId> {
        self.is_indexed = false;
        edge_table_add_row(
            &mut self.edges_,
            left.into(),
            right.into(),
            parent.into(),
            child.into(),
        )
    }

    /// Add a [``Site``] to the [``SiteTable``];
    ///
    /// # Parameters
    ///
    /// * `position`, the mutation position.
    /// * `ancestral_state`, the ancestral state at this site.
    ///
    /// # Notes
    ///
    /// If no `ancestral_state` is provided ([``None``]), then
    /// client code is assumed to have some default in mind.
    ///
    /// The [``u8``] type can be used to encode more complex
    /// state information.  Take a look at the unit tests
    /// for examples. The [bitfield](https://docs.rs/bitfield/)
    /// crate may also be useful.
    ///
    /// # Returns
    ///
    /// A [``SiteId``].
    ///
    /// The value will be [``SiteId::NULL``] in the case of overflow.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if any of the input
    /// are invalid.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// // ancestral state is a u9 equal to 3
    /// let id = tables.add_site(3, vec![3]).unwrap();
    /// assert_eq!(id, 0);
    /// // Recovering state can be a bit messy!
    /// assert_eq!(tables.site(id).ancestral_state.as_ref().unwrap(), &vec![3]);
    /// ```
    pub fn add_site<P: Into<Position>, A: Into<Option<Vec<u8>>>>(
        &mut self,
        position: P,
        ancestral_state: A,
    ) -> TablesResult<SiteId> {
        let p = position.into();
        if p >= self.length_ || p.0 < 0 {
            return Err(TablesError::InvalidPosition { found: p });
        }
        site_table_add_row(&mut self.sites_, p, ancestral_state.into())
    }

    /// Add a [``MutationRecord``] to the [``MutationTable``].
    ///
    /// # Parameters
    ///
    /// * `node`, the node where the mutation maps.
    /// * `key`, index of the mutation's metadata.
    /// * `site`, the id of the mutation's [``Site``].
    /// * `time`, the origin time of the mutation.
    /// * `derived_state`, the derived state of the variant.
    /// * `neutral`, [``true``] if the mutation affects fitness,
    ///              [``false``] otherwise.
    ///
    /// # Notes
    ///
    /// If no `derived_state` is provided ([``None``]), then
    /// client code is assumed to have some default in mind.
    ///
    /// The [``u8``] type can be used to encode more complex
    /// state information.  Take a look at the unit tests
    /// for examples. The [bitfield](https://docs.rs/bitfield/)
    /// crate may also be useful.
    ///
    /// # Returns
    ///
    /// A [``MutationId``].
    ///
    /// The value will be [``MutationId::NULL``] in the case ov
    /// overflow.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if any of the input
    /// are invalid.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// // derived state is a u9 equal to 3
    /// let id = tables.add_mutation(0, 1, 0, 0, vec![3], false).unwrap();
    /// assert_eq!(id, 0);
    /// // Recovering state can be a bit messy!
    /// assert_eq!(tables.mutation(id).derived_state.as_ref().unwrap(), &vec![3]);
    /// ```
    pub fn add_mutation<
        N: Into<NodeId>,
        S: Into<SiteId>,
        T: Into<Time>,
        K: Into<Option<usize>>,
        D: Into<Option<Vec<u8>>>,
    >(
        &mut self,
        node: N,
        key: K,
        site: S,
        time: T,
        derived_state: D,
        neutral: bool,
    ) -> TablesResult<MutationId> {
        mutation_table_add_row(
            &mut self.mutations_,
            node.into(),
            key.into(),
            site.into(),
            time.into(),
            derived_state.into(),
            neutral,
        )
    }

    /// Get genome length
    pub fn genome_length(&self) -> Position {
        self.length_
    }

    /// Return immutable reference to the [mutation table](type.MutationTable.html)
    pub fn mutations(&self) -> &[MutationRecord] {
        &self.mutations_
    }

    /// Return immutable reference to the [edge table](type.EdgeTable.html)
    pub fn edges(&self) -> &[Edge] {
        &self.edges_
    }

    /// Return number of edges
    pub fn num_edges(&self) -> usize {
        self.edges_.len()
    }

    /// Return number of nodes
    pub fn num_nodes(&self) -> usize {
        self.nodes_.len()
    }

    /// Return immutable reference to [node table](type.NodeTable.html)
    pub fn nodes(&self) -> &[Node] {
        &self.nodes_
    }

    /// Return the i-th [``Node``].
    pub fn node<N: Into<NodeId>>(&self, i: N) -> &Node {
        &self.nodes_[i.into().0 as usize]
    }

    /// Get a slice of nodes
    ///
    /// Returns `None` if the index is out of range.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// tables.add_node(0., 0).unwrap();
    /// tables.add_node(1., 1).unwrap();
    ///
    /// let n = tables.get_nodes(1).unwrap();
    /// assert_eq!(n.time, 1.into());
    /// let n = tables.get_nodes(0..2).unwrap();
    /// assert_eq!(n.len(), 2);
    /// assert_eq!(n[0].deme, 0);
    /// assert_eq!(n[1].deme, 1);
    /// ```
    pub fn get_nodes<I>(&self, index: I) -> Option<&<I as std::slice::SliceIndex<[Node]>>::Output>
    where
        I: std::slice::SliceIndex<[Node]>,
    {
        self.nodes_.get(index)
    }

    /// Return the i-th [``Edge``].
    pub fn edge<E: Into<EdgeId>>(&self, i: E) -> &Edge {
        &self.edges_[i.into().0 as usize]
    }

    /// Get a slice of edges
    ///
    /// Returns `None` if the index is out of range.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// tables.add_edge(0, 100, 0, 1).unwrap();
    /// tables.add_edge(0, 100, 0, 2).unwrap();
    ///
    /// let e = tables.get_edges(1).unwrap();
    /// assert_eq!(e.child, 2);
    /// let e = tables.get_edges(0..2).unwrap();
    /// assert_eq!(e.len(), 2);
    /// assert_eq!(e[0].child, 1);
    /// assert_eq!(e[1].child, 2);
    /// ```
    pub fn get_edges<I>(&self, index: I) -> Option<&<I as std::slice::SliceIndex<[Edge]>>::Output>
    where
        I: std::slice::SliceIndex<[Edge]>,
    {
        self.edges_.get(index)
    }

    /// Return the i-th [``Site``].
    pub fn site<S: Into<SiteId>>(&self, i: S) -> &Site {
        &self.sites_[i.into().0 as usize]
    }

    /// Get a slice of sites
    ///
    /// Returns `None` if the index is out of range.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// tables.add_site(9, None);
    /// tables.add_site(4, None);
    ///
    /// let s = tables.get_sites(1).unwrap();
    /// assert_eq!(s.position, 4);
    /// let s = tables.get_sites(0..2).unwrap();
    /// assert_eq!(s.len(), 2);
    /// assert_eq!(s[0].position, 9);
    /// assert_eq!(s[1].position, 4);
    /// ```
    pub fn get_sites<I>(&self, index: I) -> Option<&<I as std::slice::SliceIndex<[Site]>>::Output>
    where
        I: std::slice::SliceIndex<[Site]>,
    {
        self.sites_.get(index)
    }

    /// Return the i-th [``MutationRecord``].
    pub fn mutation<M: Into<MutationId>>(&self, i: M) -> &MutationRecord {
        &self.mutations_[i.into().0 as usize]
    }

    /// Get a slice of mutations
    ///
    /// Returns `None` if the index is out of range.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts_tables_trees::TableCollection::new(100).unwrap();
    /// tables.add_mutation(0, Some(113), 10, 0, None, true);
    /// tables.add_mutation(58, Some(114), 55, 0, None, true);
    ///
    /// let s = tables.get_mutations(1).unwrap();
    /// assert_eq!(s.node, 58);
    /// assert_eq!(s.site, 55);
    /// let s = tables.get_mutations(0..2).unwrap();
    /// assert_eq!(s.len(), 2);
    /// assert_eq!(s[0].site, 10);
    /// assert_eq!(s[1].site, 55);
    /// ```
    pub fn get_mutations<I>(
        &self,
        index: I,
    ) -> Option<&<I as std::slice::SliceIndex<[MutationRecord]>>::Output>
    where
        I: std::slice::SliceIndex<[MutationRecord]>,
    {
        self.mutations_.get(index)
    }

    /// Return immutable reference to [site table](type.SiteTable.html)
    pub fn sites(&self) -> &[Site] {
        &self.sites_
    }

    /// Provide an enumeration over the [node table](type.NodeTable.html)
    pub fn enumerate_nodes(&self) -> std::iter::Enumerate<std::slice::Iter<Node>> {
        self.nodes_.iter().enumerate()
    }

    /// Provide an enumeration over the [edge table](type.EdgeTable.html)
    pub fn enumerate_edges(&self) -> std::iter::Enumerate<std::slice::Iter<Edge>> {
        self.edges_.iter().enumerate()
    }

    /// Provide an enumeration over the [mutation table](type.MutationTable.html)
    pub fn enumerate_mutations(&self) -> std::iter::Enumerate<std::slice::Iter<MutationRecord>> {
        self.mutations_.iter().enumerate()
    }

    /// Provide an enumeration over the [site table](type.SiteTable.html)
    pub fn enumerate_sites(&self) -> std::iter::Enumerate<std::slice::Iter<Site>> {
        self.sites_.iter().enumerate()
    }

    /// Sort all tables for simplification.
    #[deprecated(since = "0.1.0", note = "use sort_tables instead")]
    pub fn sort_tables_for_simplification(&mut self) {
        self.sort_tables(TableSortingFlags::empty());
    }

    /// Sort all tables for simplification.
    pub fn sort_tables(&mut self, flags: TableSortingFlags) {
        if !flags.contains(TableSortingFlags::SKIP_EDGE_TABLE) {
            sort_edges(&self.nodes_, &mut self.edges_);
        }
        sort_mutation_table(&self.sites_, &mut self.mutations_);
        let mut sites: SiteTable = vec![];
        for m in self.mutations_.iter_mut() {
            record_site(&self.sites_, m, &mut sites);
        }
        std::mem::swap(&mut self.sites_, &mut sites);
    }

    /// Run a validation check on the tables.
    pub fn validate(&self, flags: TableValidationFlags) -> TablesResult<bool> {
        if flags.contains(TableValidationFlags::VALIDATE_EDGES) {
            validate_edge_table(self.genome_length(), &self.edges_, &self.nodes_)?;
        }
        if flags.contains(TableValidationFlags::VALIDATE_NODES) {
            validate_node_table(self.nodes())?;
        }
        if flags.contains(TableValidationFlags::VALIDATE_SITES) {
            validate_site_table(self.genome_length(), self.sites())?;
        }
        if flags.contains(TableValidationFlags::VALIDATE_MUTATIONS) {
            validate_mutation_table(self.mutations(), self.sites(), self.nodes())?;
        }
        Ok(true)
    }

    // SAFETY: the bounds are guaranteed by build_indexes
    fn sort_edge_output_order(edges: &[Edge], nodes: &[Node], edge_input_order: &mut Vec<usize>) {
        edge_input_order.sort_by(|a, b| {
            let ea = unsafe { edges.get_unchecked(*a) };
            let eb = unsafe { edges.get_unchecked(*b) };
            if ea.right == eb.right {
                let ta = unsafe { *nodes.get_unchecked(ea.parent.0 as usize) }.time;
                let tb = unsafe { *nodes.get_unchecked(eb.parent.0 as usize) }.time;
                match ta.partial_cmp(&tb) {
                    Some(x) => match x {
                        std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
                        std::cmp::Ordering::Less => std::cmp::Ordering::Less,
                        std::cmp::Ordering::Equal => match ea.parent.cmp(&eb.parent).reverse() {
                            std::cmp::Ordering::Equal => ea.child.cmp(&eb.child).reverse(),
                            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
                            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
                        },
                    },
                    None => panic!("invalid parent times"),
                }
            } else {
                ea.right.cmp(&eb.right)
            }
        });
    }

    // SAFETY: the bounds are guaranteed by build_indexes
    fn sort_edge_input_order(edges: &[Edge], nodes: &[Node], edge_output_order: &mut Vec<usize>) {
        edge_output_order.sort_by(|a, b| {
            let ea = unsafe { edges.get_unchecked(*a) };
            let eb = unsafe { edges.get_unchecked(*b) };
            if ea.left == eb.left {
                let ta = unsafe { *nodes.get_unchecked(ea.parent.0 as usize) }.time;
                let tb = unsafe { *nodes.get_unchecked(eb.parent.0 as usize) }.time;
                match ta.partial_cmp(&tb) {
                    Some(x) => match x.reverse() {
                        std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
                        std::cmp::Ordering::Less => std::cmp::Ordering::Less,
                        std::cmp::Ordering::Equal => match ea.parent.cmp(&eb.parent) {
                            std::cmp::Ordering::Equal => ea.child.cmp(&eb.child),
                            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
                            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
                        },
                    },
                    None => panic!("invalid parent times"),
                }
            } else {
                ea.left.cmp(&eb.left)
            }
        });
    }

    /// Build table indexes
    ///
    /// # Parameters
    ///
    /// * `flags`, see [`IndexTablesFlags`].
    ///
    /// # Errors
    ///
    /// [`TablesError`] if the input data are invalid.
    pub fn build_indexes(&mut self, flags: IndexTablesFlags) -> TablesResult<()> {
        if self.edges_.is_empty() {
            self.is_indexed = false;
            return Ok(());
        }
        if !flags.contains(IndexTablesFlags::NO_VALIDATION) {
            match validate_edge_table(self.genome_length(), &self.edges_, &self.nodes_) {
                Ok(_) => (),
                Err(e) => return Err(e),
            };
        }
        self.edge_input_order.clear();
        self.edge_output_order.clear();
        for (i, e) in self.edges_.iter().enumerate() {
            if e.parent == NodeId::NULL {
                return Err(TablesError::NullParent);
            }
            if e.child == NodeId::NULL {
                return Err(TablesError::NullChild);
            }
            if e.parent >= self.nodes_.len() as TablesIdInteger {
                return Err(TablesError::NodeOutOfBounds);
            }
            if e.child >= self.nodes_.len() as TablesIdInteger {
                return Err(TablesError::NodeOutOfBounds);
            }
            self.edge_input_order.push(i);
            self.edge_output_order.push(i);
        }
        Self::sort_edge_input_order(&self.edges_, &self.nodes_, &mut self.edge_input_order);
        Self::sort_edge_output_order(&self.edges_, &self.nodes_, &mut self.edge_output_order);
        self.is_indexed = true;
        Ok(())
    }

    /// Get the edge input order.
    ///
    /// The input order is generated by [`TableCollection::build_indexes`].
    ///
    /// Returns `None` if `self.is_indexed() == false`.
    pub fn edge_input_order(&self) -> Option<&[usize]> {
        if self.is_indexed {
            Some(&self.edge_input_order)
        } else {
            None
        }
    }

    /// Get the edge output order.
    ///
    /// The output order is generated by [`TableCollection::build_indexes`].
    ///
    /// Returns `None` if `self.is_indexed() == false`.
    pub fn edge_output_order(&self) -> Option<&[usize]> {
        if self.is_indexed {
            Some(&self.edge_output_order)
        } else {
            None
        }
    }

    /// Return `true` if tables are indexed, `false` otherwise.
    pub fn is_indexed(&self) -> bool {
        self.is_indexed
    }

    /// Dump contents of node table.
    ///
    /// The `self` object is left with an empty
    /// node table.
    #[deprecated(since = "0.3.0", note = "use dump_node_table instead")]
    pub fn dump_nodes(&mut self) -> NodeTable {
        self.dump_node_table()
    }

    /// Dump contents of edge table.
    ///
    /// The `self` object is left with an empty
    /// edge table.
    #[deprecated(since = "0.3.0", note = "use dump_edge_table instead")]
    pub fn dump_edges(&mut self) -> EdgeTable {
        self.dump_edge_table()
    }

    /// Dump contents of site table.
    ///
    /// The `self` object is left with an empty
    /// site table.
    #[deprecated(since = "0.3.0", note = "use dump_site_table instead")]
    pub fn dump_sites(&mut self) -> SiteTable {
        self.dump_site_table()
    }

    /// Dump contents of mutation table.
    ///
    /// The `self` object is left with an empty
    /// mutation table.
    #[deprecated(since = "0.3.0", note = "use dump_mutation_table instead")]
    pub fn dump_mutations(&mut self) -> MutationTable {
        self.dump_mutation_table()
    }

    /// Dump contents of node table.
    ///
    /// The `self` object is left with an empty
    /// node table.
    pub fn dump_node_table(&mut self) -> NodeTable {
        let mut rv = vec![];
        std::mem::swap(&mut self.nodes_, &mut rv);
        self.is_indexed = false;
        rv
    }

    /// Set the contents of the node table.
    pub fn set_node_table(&mut self, nodes: NodeTable) {
        self.nodes_ = nodes;
        self.is_indexed = false;
    }

    /// Dump contents of edge table.
    ///
    /// The `self` object is left with an empty
    /// edge table.
    pub fn dump_edge_table(&mut self) -> EdgeTable {
        let mut rv = vec![];
        std::mem::swap(&mut self.edges_, &mut rv);
        self.is_indexed = false;
        rv
    }

    /// Set the contents of the edge table.
    pub fn set_edge_table(&mut self, edges: EdgeTable) {
        self.edges_ = edges;
        self.is_indexed = false;
    }

    /// Dump contents of site table.
    ///
    /// The `self` object is left with an empty
    /// site table.
    pub fn dump_site_table(&mut self) -> SiteTable {
        let mut rv = vec![];
        std::mem::swap(&mut self.sites_, &mut rv);
        rv
    }

    /// Set the contents of the site table.
    pub fn set_site_table(&mut self, sites: SiteTable) {
        self.sites_ = sites;
    }

    /// Dump contents of mutation table.
    ///
    /// The `self` object is left with an empty
    /// mutation table.
    pub fn dump_mutation_table(&mut self) -> MutationTable {
        let mut rv = vec![];
        std::mem::swap(&mut self.mutations_, &mut rv);
        rv
    }

    /// Set the contents of the mutation table.
    pub fn set_mutation_table(&mut self, mutations: MutationTable) {
        self.mutations_ = mutations;
    }

    /// Count number of trees in O(E) time, where E
    /// is length of edge table.
    ///
    /// # Errors
    ///
    /// [`TablesError::TablesNotIndexed`] if tables are not indexed
    ///
    /// # Panics
    ///
    /// If the edge table is invalid in any way, a `panic!` may occur.
    /// To check table validity, call [`TableCollection::validate`].
    pub fn count_trees(&self) -> TablesResult<u32> {
        if !self.is_indexed() {
            Err(TablesError::TablesNotIndexed)
        } else {
            let mut num_trees = 0;
            let mut input_index: usize = 0;
            let mut output_index: usize = 0;
            let input = self.edge_input_order.as_slice();
            let output = self.edge_output_order.as_slice();
            let edges = self.edges_.as_slice();

            let mut tree_left = Position(0);
            while input_index < input.len() || tree_left < self.genome_length() {
                for idx in output[output_index..].iter() {
                    if edges[*idx].right != tree_left {
                        break;
                    }
                    output_index += 1;
                }
                for idx in input[input_index..].iter() {
                    if edges[*idx].left != tree_left {
                        break;
                    }
                    input_index += 1;
                }
                let mut tree_right = self.genome_length();
                if input_index < input.len() {
                    tree_right = std::cmp::min(tree_right, edges[input[input_index]].left);
                }
                if output_index < output.len() {
                    tree_right = std::cmp::min(tree_right, edges[output[output_index]].right);
                }
                tree_left = tree_right;
                num_trees += 1;
            }
            Ok(num_trees)
        }
    }
}

#[cfg(test)]
mod test_tables {

    use super::*;

    #[test]
    fn test_flags_vs_tskit() {
        assert_eq!(NodeFlags::IS_SAMPLE.bits(), tskit::TSK_NODE_IS_SAMPLE);
    }

    #[test]
    fn test_bad_genome_length() {
        let _ = TableCollection::new(Position(0)).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidGenomeLength),
            |_| panic!(),
        );
    }

    #[test]
    fn test_add_edge() {
        let mut tables = TableCollection::new(10).unwrap();

        let result = tables.add_edge(0, 1, 2, 3).unwrap();

        assert_eq!(0, result);
        assert_eq!(1, tables.edges().len());
        assert_eq!(1, tables.num_edges());

        let mut n = 0;
        for _ in tables.edges() {
            n += 1;
        }
        assert_eq!(n, 1);
    }

    #[test]
    fn test_add_edge_bad_positions() {
        let mut tables = TableCollection::new(10).unwrap();

        let _ = tables.add_edge(-1, 1, 1, 2).map_or_else(
            |x: TablesError| {
                assert_eq!(
                    x,
                    TablesError::InvalidPosition {
                        found: Position(-1)
                    }
                )
            },
            |_| panic!(),
        );

        let _ = tables.add_edge(1, -1, 1, 2).map_or_else(
            |x: TablesError| {
                assert_eq!(
                    x,
                    TablesError::InvalidLeftRight {
                        found: (Position(1), Position(-1))
                    }
                )
            },
            |_| panic!(),
        );
    }

    #[test]
    fn test_add_edge_bad_nodes() {
        let mut tables = TableCollection::new(10).unwrap();

        let _ = tables.add_edge(0, 1, -1, 2).map_or_else(
            |x: TablesError| {
                assert_eq!(
                    x,
                    TablesError::InvalidNodeValue {
                        found: NodeId::NULL
                    }
                )
            },
            |_| panic!(),
        );

        let _ = tables.add_edge(0, 1, 1, -2).map_or_else(
            |x: TablesError| {
                assert_eq!(
                    x,
                    TablesError::InvalidNodeValue {
                        found: NodeId::NULL
                    }
                )
            },
            |_| panic!(),
        );
    }

    #[test]
    #[should_panic]
    fn test_add_site_negative_position() {
        let mut tables = TableCollection::new(10).unwrap();
        tables.add_site(-1, None).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_add_site_position_too_big() {
        let mut tables = TableCollection::new(10).unwrap();
        tables.add_site(tables.genome_length(), None).unwrap();
    }

    #[test]
    fn test_add_site_with_ancestral_state() {
        let mut tables = TableCollection::new(10).unwrap();
        let _ = tables
            .add_site(1, Some(b"0".to_vec()))
            .map_or_else(|_: TablesError| panic!(), |_| ());
        let s = tables.site(0);
        assert_eq!(s.position, 1);
        assert_eq!(s.ancestral_state, Some(b"0".to_vec()));
        assert_eq!(s.ancestral_state.as_ref().unwrap().len(), 1);

        match std::str::from_utf8(s.ancestral_state.as_ref().unwrap()) {
            Ok(x) => assert_eq!(x, "0".to_string()),
            Err(_) => panic!(),
        }
    }

    fn decode_complex_state<T>(v: &[T]) -> &[T; 4] {
        use std::convert::TryInto;
        v.try_into().unwrap_or_else(|_| panic!())
    }

    #[test]
    fn test_add_site_with_complex_ancestral_state() {
        let mut tables = TableCollection::new(Position(10)).unwrap();
        let mut a: Vec<u8> = vec![];
        let astate_part1: i32 = 300;
        let astate_part2: i32 = 3;
        a.append(&mut astate_part1.to_le_bytes().to_vec());
        a.append(&mut astate_part2.to_le_bytes().to_vec());
        assert_eq!(a.len(), 8);
        let _ = tables
            .add_site(1, Some(a))
            .map_or_else(|_: TablesError| panic!(), |_| ());

        // now, decode it
        let s = tables.site(0);

        // The conversion works by dereferencing a reference
        // to a mutable array whose lifetime is tied to the
        // slice:
        let x = i32::from_le_bytes(*decode_complex_state(
            &s.ancestral_state.as_ref().unwrap()[0..4],
        ));
        assert_eq!(x, astate_part1);
        let x = i32::from_le_bytes(*decode_complex_state(
            &s.ancestral_state.as_ref().unwrap()[4..8],
        ));
        assert_eq!(x, astate_part2);
    }

    #[test]
    fn test_add_site_without_ancestral_state() {
        let mut tables = TableCollection::new(10).unwrap();
        let _ = tables
            .add_site(1, None)
            .map_or_else(|_: TablesError| panic!(), |_| ());
        let s = tables.site(0);
        if s.ancestral_state.as_ref().is_some() {
            panic!()
        }
    }

    #[test]
    fn test_add_mutation_without_derived_state() {
        let mut tables = TableCollection::new(10).unwrap();
        let _ = tables.add_mutation(0, None, 0, 0, None, false).unwrap();
        let m = tables.mutation(0);
        if m.derived_state.as_ref().is_some() {
            panic!()
        }
    }

    #[test]
    fn test_add_mutation_with_derived_state() {
        let mut tables = TableCollection::new(10).unwrap();
        let _ = tables
            .add_mutation(0, None, 0, 0, Some(b"0".to_vec()), false)
            .unwrap();
        let m = tables.mutation(0);
        match std::str::from_utf8(m.derived_state.as_ref().unwrap()) {
            Ok(x) => assert_eq!(x, "0".to_string()),
            Err(_) => panic!(),
        }
    }

    #[test]
    #[allow(clippy::redundant_clone)]
    fn test_clone_tables() {
        let mut tables = TableCollection::new(10).unwrap();
        tables.add_edge(0, 5, 0, 1).unwrap();
        let tclone = tables.clone();

        assert_eq!(tclone.edges().len(), 1);
        let e = tclone.edge(0);
        assert_eq!(e.left, 0);
        assert_eq!(e.right, 5);
        assert_eq!(e.parent, 0);
        assert_eq!(e.child, 1);
    }

    #[test]
    fn test_node_flags() {
        let mut x = (NodeFlags::IS_ALIVE | NodeFlags::IS_SAMPLE).bits();
        assert!(x & NodeFlags::IS_ALIVE.bits() > 0);
        assert!(x & NodeFlags::IS_SAMPLE.bits() > 0);
        x &= !NodeFlags::IS_SAMPLE.bits();
        assert!(x & NodeFlags::IS_ALIVE.bits() > 0);
        assert!(x & NodeFlags::IS_SAMPLE.bits() == 0);
    }
}

#[cfg(test)]
mod test_table_indexing {
    use super::*;

    #[test]
    fn test_reverse_sort() {
        let mut v = vec![3, 2, 1, 0];
        v.sort_by(|a, b| a.cmp(b).reverse());
        for i in 1..v.len() {
            assert!(v[i] <= v[i - 1]);
        }
    }

    #[test]
    #[should_panic]
    fn test_no_nodes() {
        let mut t = TableCollection::new(1).unwrap();
        t.add_edge(0, 1, 0, 1).unwrap();
        t.add_edge(0, 1, 0, 2).unwrap();
        t.build_indexes(IndexTablesFlags::default()).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_edge_out_of_range() {
        let mut t = TableCollection::new(1).unwrap();
        t.add_node(0., 0).unwrap();
        t.add_edge(0, 1, 0, 1).unwrap();
        t.build_indexes(IndexTablesFlags::default()).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_simple_invalid_edge_table() {
        let mut t = TableCollection::new(1).unwrap();
        for _ in 0..3 {
            t.add_node(2., 0).unwrap();
        }
        t.add_node(1., 0).unwrap();
        t.add_node(0., 0).unwrap();

        t.add_edge(0, 1, 4, 3).unwrap();
        t.add_edge(0, 1, 4, 2).unwrap();
        t.add_edge(0, 1, 3, 0).unwrap();
        t.add_edge(0, 1, 3, 1).unwrap();

        assert_eq!(t.nodes().len(), 5);

        validate_edge_table(t.genome_length(), t.edges(), t.nodes()).unwrap();
        t.build_indexes(IndexTablesFlags::default()).unwrap();
    }

    #[test]
    fn test_simple_sort_order() {
        let mut t = TableCollection::new(1).unwrap();
        for _ in 0..3 {
            t.add_node(2., 0).unwrap();
        }
        t.add_node(1., 0).unwrap();
        t.add_node(0., 0).unwrap();

        t.add_edge(0, 1, 4, 3).unwrap();
        t.add_edge(0, 1, 4, 2).unwrap();
        t.add_edge(0, 1, 3, 0).unwrap();
        t.add_edge(0, 1, 3, 1).unwrap();

        assert_eq!(t.nodes().len(), 5);

        t.sort_tables(TableSortingFlags::empty());
        validate_edge_table(t.genome_length(), t.edges(), t.nodes()).unwrap();
        t.build_indexes(IndexTablesFlags::default()).unwrap();

        if let Some(edge_input_order) = t.edge_input_order() {
            assert_eq!(edge_input_order.len(), t.edges().len());
            for (idx, i) in edge_input_order.iter().enumerate() {
                if idx > 0 {
                    let ti = t.node(t.edge(EdgeId::from(*i)).parent).time;
                    let tim1 = t
                        .node(t.edge(EdgeId::from(edge_input_order[idx] - 1)).parent)
                        .time;
                    assert!(ti <= tim1);
                    assert!(tim1 >= ti);
                }
            }
        } else {
            panic!("expected a edge_input_order");
        }

        if let Some(edge_output_order) = t.edge_output_order() {
            assert_eq!(edge_output_order.len(), t.edges().len());
            for (idx, i) in edge_output_order.iter().enumerate() {
                if idx > 0 {
                    let ti = t.node(t.edge(EdgeId::from(*i)).parent).time;
                    let tim1 = t
                        .node(t.edge(EdgeId::from(edge_output_order[idx - 1])).parent)
                        .time;
                    assert!(ti >= tim1, "{} {}", f64::from(ti), f64::from(tim1));
                }
            }
        } else {
            panic!("expected a edge_output_order");
        }
    }

    #[test]
    fn test_is_indexed() {
        let mut t = TableCollection::new(1).unwrap();
        for _ in 0..3 {
            t.add_node(2., 0).unwrap();
        }
        t.add_node(1., 0).unwrap();
        t.add_node(0., 0).unwrap();

        t.add_edge(0, 1, 4, 3).unwrap();
        t.add_edge(0, 1, 4, 2).unwrap();
        t.add_edge(0, 1, 3, 0).unwrap();
        t.add_edge(0, 1, 3, 1).unwrap();

        t.sort_tables(TableSortingFlags::empty());
        validate_edge_table(t.genome_length(), t.edges(), t.nodes()).unwrap();
        t.build_indexes(IndexTablesFlags::default()).unwrap();

        assert!(t.is_indexed());

        t.add_edge(0, 1, 4, 0).unwrap();

        assert!(!t.is_indexed());

        t.sort_tables(TableSortingFlags::empty());
        validate_edge_table(t.genome_length(), t.edges(), t.nodes()).unwrap();
        t.build_indexes(IndexTablesFlags::default()).unwrap();
        assert!(t.is_indexed());

        t.add_node(0., 0).unwrap();
        assert!(!t.is_indexed());
    }
}

#[cfg(test)]
mod test_table_validation {
    use super::*;

    #[test]
    fn test_validation_flags() {
        let v = vec![
            TableValidationFlags::VALIDATE_EDGES,
            TableValidationFlags::VALIDATE_SITES,
            TableValidationFlags::VALIDATE_MUTATIONS,
        ];
        for f in v.iter() {
            for ff in v.iter() {
                if *f != *ff {
                    assert!(!f.contains(*ff));
                }
            }
        }
    }

    #[test]
    fn test_site_table_not_sorted_by_position() {
        // edges aren't sorted, but we skip that check
        let mut t = TableCollection::new(10).unwrap();
        let node0 = t.add_node(0., 0).unwrap();
        let node1 = t.add_node(1., 0).unwrap();
        t.add_edge(0, t.genome_length(), node1, node0).unwrap();
        t.add_site(5, None).unwrap();
        t.add_site(4, None).unwrap();
        match t.validate(TableValidationFlags::VALIDATE_SITES) {
            Err(TablesError::UnsortedSitePosition) => (),
            Err(_) => panic!("unexpected Err"),
            Ok(_) => panic!("unexpected Ok"),
        };
    }
}
