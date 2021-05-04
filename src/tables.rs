use crate::tsdef::{IdType, Position, Time, NULL_ID};
use bitflags::bitflags;
use std::cmp::Ordering;
use thiserror::Error;

/// Error type related to [``TableCollection``]
#[derive(Error, Debug, PartialEq)]
pub enum TablesError {
    /// Raised by [``TableCollection::new``].
    #[error("Invalid genome length")]
    InvalidGenomeLength,
    /// Raised when invalid node `ID`s are encountered.
    #[error("Invalid node: {found:?}")]
    InvalidNodeValue {
        /// The invalid `ID`
        found: IdType,
    },
    /// Raised when invalid positions are encountered.
    #[error("Invalid value for position: {found:?}")]
    InvalidPosition {
        /// The invalid position
        found: Position,
    },
    /// Raised when an [``Edge``]'s left/right
    /// values are invalid.
    #[error("Invalid position range: {found:?}")]
    InvalidLeftRight {
        /// The invalid `(left, right)`.
        found: (Position, Position),
    },
    /// Raised when invalid times are encountered.
    #[error("Invalid value for time: {found:?}")]
    InvalidTime {
        /// The invalid time
        found: Time,
    },
    #[error("Invalid value for deme: {found:?}")]
    /// Raised with a deme's `ID` is invalid.
    InvalidDeme {
        /// The invalide deme `ID`
        found: IdType,
    },
    #[error("Parent is NULL_ID")]
    /// Can be raised by [``validate_edge_table``]
    NullParent,
    #[error("Child is NULL_ID")]
    /// Can be raised by [``validate_edge_table``]
    NullChild,
    #[error("Node is out of bounds")]
    /// Can be raised by [``validate_edge_table``]
    NodeOutOfBounds,
    #[error("Node time order violation")]
    /// Can be raised by [``validate_edge_table``]
    NodeTimesUnordered,
    #[error("Parents not sorted by time")]
    /// Can be raised by [``validate_edge_table``]
    ParentTimesUnsorted,
    #[error("Parents not contiguous")]
    /// Can be raised by [``validate_edge_table``]
    ParentsNotContiguous,
    #[error("Edges not sorted by child")]
    /// Can be raised by [``validate_edge_table``]
    EdgesNotSortedByChild,
    #[error("Edges not sorted by left")]
    /// Can be raised by [``validate_edge_table``]
    EdgesNotSortedByLeft,
    #[error("Duplicate edges")]
    /// Can be raised by [``validate_edge_table``]
    DuplicateEdges,
}

/// Result type for operations on tables
pub type TablesResult<T> = std::result::Result<T, TablesError>;

/// A Node of a tree sequence
#[derive(Copy, Clone)]
pub struct Node {
    /// Birth time
    pub time: Time,
    /// Population (deme) of node
    pub deme: IdType,
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
    pub parent: IdType,
    /// Index of child in a [NodeTable](type.NodeTable.html)
    pub child: IdType,
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
    pub node: IdType,
    /// Reference to mutation metadata.
    pub key: usize,
    /// The index of the corresponding [``Site``].
    pub site: usize,
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
    if x < 0 {
        Err(TablesError::InvalidPosition { found: x })
    } else {
        Ok(())
    }
}

fn node_non_negative(x: IdType) -> TablesResult<()> {
    if x < 0 {
        Err(TablesError::InvalidNodeValue { found: x })
    } else {
        Ok(())
    }
}

//fn time_non_negative(x: Time) -> TablesResult<()> {
//    if x < 0 {
//        Err(TablesError::InvalidTime { found: x })
//    } else {
//        Ok(())
//    }
//}

fn deme_non_negative(x: IdType) -> TablesResult<()> {
    if x < 0 {
        Err(TablesError::InvalidDeme { found: x })
    } else {
        Ok(())
    }
}

fn edge_table_add_row(
    edges: &mut EdgeTable,
    left: Position,
    right: Position,
    parent: IdType,
    child: IdType,
) -> TablesResult<IdType> {
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

    Ok((edges.len() - 1) as IdType)
}

// NOTE: we allow negative times, in order to support "precapitation".
fn node_table_add_row(nodes: &mut NodeTable, time: Time, deme: IdType) -> TablesResult<IdType> {
    //time_non_negative(time)?;
    deme_non_negative(deme)?;
    nodes.push(Node { time, deme });

    Ok((nodes.len() - 1) as IdType)
}

fn site_table_add_row(
    sites: &mut SiteTable,
    position: Position,
    ancestral_state: Option<Vec<u8>>,
) -> TablesResult<IdType> {
    position_non_negative(position)?;
    sites.push(Site {
        position,
        ancestral_state,
    });
    Ok((sites.len() - 1) as IdType)
}

fn mutation_table_add_row(
    mutations: &mut MutationTable,
    node: IdType,
    key: usize,
    site: usize,
    derived_state: Option<Vec<u8>>,
    neutral: bool,
) -> TablesResult<IdType> {
    node_non_negative(node)?;
    mutations.push(MutationRecord {
        node,
        key,
        site,
        derived_state,
        neutral,
    });
    Ok((mutations.len() - 1) as IdType)
}

fn sort_edges(nodes: &[Node], edges: &mut [Edge]) {
    edges.sort_by(|a, b| {
        let aindex = a.parent as usize;
        let bindex = b.parent as usize;
        let ta = nodes[aindex].time;
        let tb = nodes[bindex].time;
        if ta == tb {
            if a.parent == b.parent {
                if a.child == b.child {
                    return a.left.cmp(&b.left);
                }
                return a.child.cmp(&b.child);
            }
            return a.parent.cmp(&b.parent);
        }
        ta.cmp(&tb).reverse()
    });
}

fn record_site(sites: &[Site], mutation: &mut MutationRecord, new_site_table: &mut SiteTable) {
    let position = sites[mutation.site].position;
    if new_site_table.is_empty() || new_site_table[new_site_table.len() - 1].position != position {
        new_site_table.push(sites[mutation.site as usize].clone());
    }

    mutation.site = new_site_table.len() - 1;
}

fn sort_mutation_table(sites: &[Site], mutations: &mut [MutationRecord]) {
    mutations.sort_by(|a, b| {
        let pa = sites[a.site].position;
        let pb = sites[b.site].position;
        pa.cmp(&pb)
    });
}

bitflags! {
    /// Modifies behavior of
    /// [``TableCollection::validate``]
    ///
    /// ```
    /// let f = forrustts::TableValidationFlags::empty();
    /// assert_eq!(f.contains(forrustts::TableValidationFlags::VALIDATE_ALL), true);
    /// ```
    #[derive(Default)]
    pub struct TableValidationFlags: u32 {
        /// Validate all tables.
        /// This is also the "default"/empty.
        const VALIDATE_ALL = 0;
    }
}

bitflags! {
    /// Modifies behavior of
    /// [``TableCollection::sort_tables``]
    ///
    /// ```
    /// let f = forrustts::TableSortingFlags::empty();
    /// assert_eq!(f.contains(forrustts::TableSortingFlags::SORT_ALL), true);
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
/// let mut tables = forrustts::TableCollection::new(100).unwrap();
/// // (do some stuff now...)
/// let rv = forrustts::validate_edge_table(tables.genome_length(),
///                                         &tables.edges(),
///                                         &tables.nodes()).unwrap();
/// assert_eq!(rv, true);
/// ```
pub fn validate_edge_table(len: Position, edges: &[Edge], nodes: &[Node]) -> TablesResult<bool> {
    if edges.is_empty() {
        return Ok(true);
    }
    let mut parent_seen = vec![0; nodes.len()];
    let mut last_parent: usize = edges[0].parent as usize;
    let mut last_child: usize = edges[0].child as usize;
    let mut last_left: Position = edges[0].left;

    for (i, edge) in edges.iter().enumerate() {
        if edge.parent == NULL_ID {
            return Err(TablesError::NullParent);
        }
        if edge.child == NULL_ID {
            return Err(TablesError::NullChild);
        }
        if edge.parent < 0 || edge.parent as usize >= nodes.len() {
            return Err(TablesError::NodeOutOfBounds);
        }
        if edge.child < 0 || edge.child as usize >= nodes.len() {
            return Err(TablesError::NodeOutOfBounds);
        }
        if edge.left < 0 || edge.left > len {
            return Err(TablesError::InvalidPosition { found: edge.left });
        }
        if edge.right < 0 || edge.right > len {
            return Err(TablesError::InvalidPosition { found: edge.right });
        }
        if edge.left >= edge.right {
            return Err(TablesError::InvalidLeftRight {
                found: (edge.left, edge.right),
            });
        }

        // child time must be > parent time b/c time goes forwards
        if nodes[edge.child as usize].time <= nodes[edge.parent as usize].time {
            return Err(TablesError::NodeTimesUnordered);
        }

        if parent_seen[edge.parent as usize] == 1 {
            return Err(TablesError::ParentsNotContiguous);
        }

        if i > 0 {
            if nodes[edge.parent as usize].time > nodes[last_parent].time {
                return Err(TablesError::ParentTimesUnsorted);
            }
            if nodes[edge.parent as usize].time == nodes[last_parent].time {
                if edge.parent as usize == last_parent {
                    if (edge.child as usize) < last_child {
                        return Err(TablesError::EdgesNotSortedByChild);
                    }
                    if edge.child as usize == last_child {
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
        }
        last_parent = edge.parent as usize;
        last_child = edge.child as usize;
        last_left = edge.left;
    }

    Ok(true)
}

/// A collection of node, edge, site, and mutation tables.
#[derive(Clone)]
pub struct TableCollection {
    length_: Position, // Not visible outside of this module

    pub(crate) nodes_: NodeTable,
    pub(crate) edges_: EdgeTable,
    pub(crate) sites_: SiteTable,
    pub(crate) mutations_: MutationTable,
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
    pub const fn new(genome_length: Position) -> TablesResult<TableCollection> {
        if genome_length < 1 {
            return Err(TablesError::InvalidGenomeLength);
        }

        Ok(TableCollection {
            length_: genome_length,
            nodes_: NodeTable::new(),
            edges_: EdgeTable::new(),
            sites_: SiteTable::new(),
            mutations_: MutationTable::new(),
        })
    }

    /// Add a [``Node``] to the [``NodeTable``]
    ///
    /// # Parameters
    ///
    /// * `time`, a [``Time``] representing the birth time.
    /// * `deme` a valid [``IdType``] representing deme where the node is found.
    ///
    /// # Returns
    ///
    /// An [``IdType``] that is the new node's ``ID``.
    /// This value is the index of the node in the node table.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if `deme < 0`.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts::TableCollection::new(100).unwrap();
    /// let id = tables.add_node(1, 0).unwrap();
    /// assert_eq!(id, 0);
    /// ```
    pub fn add_node(&mut self, time: Time, deme: IdType) -> TablesResult<IdType> {
        node_table_add_row(&mut self.nodes_, time, deme)
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
    /// An [``IdType``] that is the new edge's ``ID``.
    /// This value is the index of the edge in the edge table.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if any of the input
    /// are invalid.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts::TableCollection::new(100).unwrap();
    /// let id = tables.add_edge(0, 3, 5, 9).unwrap();
    /// assert_eq!(id, 0);
    /// ```
    pub fn add_edge(
        &mut self,
        left: Position,
        right: Position,
        parent: IdType,
        child: IdType,
    ) -> TablesResult<IdType> {
        edge_table_add_row(&mut self.edges_, left, right, parent, child)
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
    /// An [``IdType``] that is the new site's ``ID``.
    /// This value is the index of the site in the site table.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if any of the input
    /// are invalid.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts::TableCollection::new(100).unwrap();
    /// // ancestral state is a u9 equal to 3
    /// let id = tables.add_site(3, Some(vec![3])).unwrap();
    /// assert_eq!(id, 0);
    /// // Recovering state can be a bit messy!
    /// assert_eq!(tables.site(id).ancestral_state.as_ref().unwrap(), &vec![3]);
    /// ```
    pub fn add_site(
        &mut self,
        position: Position,
        ancestral_state: Option<Vec<u8>>,
    ) -> TablesResult<IdType> {
        if position >= self.length_ || position < 0 {
            return Err(TablesError::InvalidPosition { found: position });
        }
        site_table_add_row(&mut self.sites_, position, ancestral_state)
    }

    /// Add a [``MutationRecord``] to the [``MutationTable``].
    ///
    /// # Parameters
    ///
    /// * `node`, the node where the mutation maps.
    /// * `key`, index of the mutation's metadata.
    /// * `site`, the [``IdType``] of the mutation's [``Site``].
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
    /// An [``IdType``] that is the new mutation's ``ID``.
    /// This value is the index of the mutation in the mutation table.
    ///
    /// # Errors
    ///
    /// Will return [``TablesError``] if any of the input
    /// are invalid.
    ///
    /// # Example
    ///
    /// ```
    /// let mut tables = forrustts::TableCollection::new(100).unwrap();
    /// // derived state is a u9 equal to 3
    /// let id = tables.add_mutation(0, 0, 0, Some(vec![3]), false).unwrap();
    /// assert_eq!(id, 0);
    /// // Recovering state can be a bit messy!
    /// assert_eq!(tables.mutation(id).derived_state.as_ref().unwrap(), &vec![3]);
    /// ```
    pub fn add_mutation(
        &mut self,
        node: IdType,
        key: usize,
        site: usize,
        derived_state: Option<Vec<u8>>,
        neutral: bool,
    ) -> TablesResult<IdType> {
        mutation_table_add_row(
            &mut self.mutations_,
            node,
            key,
            site,
            derived_state,
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
    pub fn node(&self, i: IdType) -> &Node {
        &self.nodes_[i as usize]
    }

    /// Return the i-th [``Edge``].
    pub fn edge(&self, i: IdType) -> &Edge {
        &self.edges_[i as usize]
    }

    /// Return the i-th [``Site``].
    pub fn site(&self, i: IdType) -> &Site {
        &self.sites_[i as usize]
    }

    /// Return the i-th [``MutationRecord``].
    pub fn mutation(&self, i: IdType) -> &MutationRecord {
        &self.mutations_[i as usize]
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
        if flags.contains(TableValidationFlags::VALIDATE_ALL) {
            validate_edge_table(self.genome_length(), &self.edges_, &self.nodes_)?;
        }
        Ok(true)
    }
}

#[cfg(test)]
mod test_tables {

    use super::*;

    #[test]
    fn test_bad_genome_length() {
        let _ = TableCollection::new(0).map_or_else(
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
            |x: TablesError| assert_eq!(x, TablesError::InvalidPosition { found: -1 }),
            |_| panic!(),
        );

        let _ = tables.add_edge(1, -1, 1, 2).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidLeftRight { found: (1, -1) }),
            |_| panic!(),
        );
    }

    #[test]
    fn test_add_edge_bad_nodes() {
        let mut tables = TableCollection::new(10).unwrap();

        let _ = tables.add_edge(0, 1, -1, 2).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidNodeValue { found: -1 }),
            |_| panic!(),
        );

        let _ = tables.add_edge(0, 1, 1, -2).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidNodeValue { found: -2 }),
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
        let mut tables = TableCollection::new(10).unwrap();
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
        let _ = tables.add_mutation(0, 0, 0, None, false).unwrap();
        let m = tables.mutation(0);
        if m.derived_state.as_ref().is_some() {
            panic!()
        }
    }

    #[test]
    fn test_add_mutation_with_derived_state() {
        let mut tables = TableCollection::new(10).unwrap();
        let _ = tables
            .add_mutation(0, 0, 0, Some(b"0".to_vec()), false)
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
}
