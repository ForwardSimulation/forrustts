use crate::tsdef::{POSITION, TIME};
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum TablesError {
    #[error("Invalid genome length")]
    InvalidGenomeLength,
    #[error("Invalid node: {found:?}")]
    InvalidNodeValue { found: i64 }, // FIXME: is there a way to make this generic?
    #[error("Invalid value for position: {found:?}")]
    InvalidPosition { found: i64 },
    #[error("Invalid position range: {found:?}")]
    InvalidLeftRight { found: (i64, i64) },
    #[error("Invalid value for time: {found:?}")]
    InvalidTime { found: i64 },
    #[error("Invalid value for deme: {found:?}")]
    InvalidDeme { found: i32 },
    #[error("Parent is NULLTSINT")]
    NullParent,
    #[error("Child is NULLTSINT")]
    NullChild,
    #[error("Node is out of bounds")]
    NodeOutOfBounds,
    #[error("Node time order violation")]
    NodeTimesUnordered,
    #[error("Parents not sorted by time")]
    ParentTimesUnsorted,
    #[error("Parents not contiguous")]
    ParentsNotContiguous,
    #[error("Edges not sorted by child")]
    EdgesNotSortedByChild,
    #[error("Edges not sorted by left")]
    EdgesNotSortedByLeft,
    #[error("Duplicate edges")]
    DuplicateEdges,
}

/// Result type for operations on tables
pub type TablesResult<T> = std::result::Result<T, TablesError>;

/// A Node of a tree sequence
pub struct Node<T> {
    /// Birth time
    pub time: TIME,
    /// Population (deme) of node
    pub deme: T,
}

/// An Edge is a transmission event
#[derive(Copy, Clone)]
pub struct Edge<T> {
    pub left: POSITION,
    pub right: POSITION,
    /// Index of parent in a [NodeTable](type.NodeTable.html)
    pub parent: T,
    /// Index of child in a [NodeTable](type.NodeTable.html)
    pub child: T,
}

// TODO: It would be nice to use generics here
// to allow arbitrary types for ancestral_state
// and derived_state.

/// A Site is the location and
/// ancestral state of a tables::Mutation
pub struct Site {
    pub position: POSITION,
    pub ancestral_state: i8,
}

/// A Mutation is the minimal information
/// needed about a mutation to track it
/// on a tree sequence.
pub struct Mutation<T> {
    pub node: T,
    pub key: usize,
    pub site: usize,
    pub derived_state: i8,
    pub neutral: bool,
}

// TODO: do these need to be pub
pub type NodeTable<T> = Vec<Node<T>>;
pub type EdgeTable<T> = Vec<Edge<T>>;
pub type SiteTable = Vec<Site>;
pub type MutationTable<T> = Vec<Mutation<T>>;

fn position_non_negative(x: POSITION) -> TablesResult<()> {
    if x < 0 {
        return Err(TablesError::InvalidPosition { found: x });
    }
    Ok(())
}

fn node_non_negative<T>(x: T) -> TablesResult<()> {
    if x < 0 {
        return Err(TablesError::InvalidNodeValue { found: x });
    }
    Ok(())
}

fn time_non_negative(x: TIME) -> TablesResult<()> {
    if x < 0 {
        return Err(TablesError::InvalidTime { found: x });
    }
    Ok(())
}

fn deme_non_negative<T>(x: T) -> TablesResult<()> {
    if x < 0 {
        return Err(TablesError::InvalidDeme { found: x });
    }
    Ok(())
}

macro_rules! add_edge {
    ($itype: ty) => {
        fn add_edge(
            edges: &mut EdgeTable<$itype>,
            left: POSITION,
            right: POSITION,
            parent: $itype,
            child: $itype,
        ) -> TablesResult<$itype> {
            if right <= left {
                return Err(TablesError::InvalidLeftRight {
                    found: (left, right),
                });
            }
            position_non_negative(left)?;
            position_non_negative(right)?;
            node_non_negative(parent)?;
            node_non_negative(child)?;

            edges.push(Edge::<$itype> {
                left: left,
                right: right,
                parent: parent,
                child: child,
            });

            Ok(edges.len() as $itype)
        }
    };
}

macro_rules! add_node {
    ($itype: ty) => {
        fn add_node(
            nodes: &mut NodeTable<$itype>,
            time: TIME,
            deme: $itype,
        ) -> TablesResult<$itype> {
            time_non_negative(time)?;
            deme_non_negative(deme)?;
            nodes.push(Node::<$itype> {
                time: time,
                deme: deme,
            });

            // TODO: learn if there is a way to raise error
            // automagically if overlow.
            Ok((nodes.len() - 1) as $itype)
        }
    };
}

macro_rules! add_site {
    ($itype: ty) => {
        fn add_site(
            sites: &mut SiteTable,
            position: POSITION,
            ancestral_state: i8,
        ) -> TablesResult<$itype> {
            position_non_negative(position)?;
            sites.push(Site {
                position: position,
                ancestral_state: ancestral_state,
            });
            Ok(sites.len() as $itype)
        }
    };
}

macro_rules! add_mutation {
    ($itype: ty) => {
        fn add_mutation(
            mutations: &mut MutationTable<$itype>,
            node: $itype,
            key: usize,
            site: usize,
            derived_state: i8,
            neutral: bool,
        ) -> TablesResult<$itype> {
            node_non_negative(node)?;
            mutations.push(Mutation::<$itype> {
                node: node,
                key: key,
                site: site,
                derived_state: derived_state,
                neutral: neutral,
            });
            Ok(mutations.len() as $itype)
        }
    };
}

fn sort_edge_table<T>(nodes: &NodeTable<T>, edges: &mut EdgeTable<T>) -> () {
    // NOTE: it may by more idiomatic to
    // not use a slice here, and instead allow
    // the range-checking?
    let nslice = &nodes.as_slice();
    edges.sort_by(|a, b| {
        let aindex = a.parent as usize;
        let bindex = b.parent as usize;
        let ta = nslice[aindex].time;
        let tb = nslice[bindex].time;
        if ta == tb {
            if a.parent == b.parent {
                if a.child == b.child {
                    return a.left.cmp(&b.left);
                }
                return a.child.cmp(&b.child);
            }
            return a.parent.cmp(&b.parent);
        }
        return ta.cmp(&tb).reverse();
    });
}

fn sort_mutation_table<T>(sites: &SiteTable, mutations: &mut MutationTable<T>) -> () {
    let sslice = &sites.as_slice();
    mutations.sort_by(|a, b| {
        let pa = sslice[a.site].position;
        let pb = sslice[b.site].position;
        return pa.partial_cmp(&pb).unwrap().reverse();
    });
}

macro_rules! sort_tables_for_simplification {
    ($itype: ty) => {
        fn sort_tables_for_simplification(&mut self) -> () {
            sort_edge_table(&self.nodes_, &mut self.edges_);
            sort_mutation_table(&self.nodes_, &mut self.edges_);
        }
    };
}

macro_rules! tree_sequence_recording_interface {
    ($itype: ty) => {
        // The public interfact
        add_node!($itype);
        add_edge!($itype);
        add_site!($itype);
        add_mutation!($itype);
        sort_tables_for_simplification!($itype);
    };
}

//FIXME: remove later.
//This is a mock-up

// Table row
pub struct E<T> {
    left: i64,
    right: i64,
    parent: T,
    child: T,
}

// table collection
pub struct TC<T> {
    edges_: Vec<E<T>>,
}

// Interface
trait TCT<T> {
    fn add_edge(&mut self, left: i64, right: i64, parent: T, child: T) -> TablesResult<usize>;
    fn foo();
}

// macros for the details
macro_rules! add_edge2 {
    ($itype: ty) => {
        fn add_edge(
            &mut self,
            left: i64,
            right: i64,
            parent: $itype,
            child: $itype,
        ) -> TablesResult<$itype> {
            self.edges_.push(E::<$itype> {
                left,
                right,
                parent,
                child,
            });
            Ok(self.edges_.len() as $itype)
        }
    };
}

macro_rules! add_foo {
    ($itype: ty) => {
        fn foo() {}
    };
}

macro_rules! mock_make_interface {
    ($itype: ty) => {
        add_edge2!($itype);
        add_foo!($itype);
    };
}

// Concrete definitions
impl TCT<i32> for TC<i32> {
    mock_make_interface!(i32);
}

impl TCT<i64> for TC<i64> {
    mock_make_interface!(i64);
}

pub fn validate_edge_table<T>(
    len: POSITION,
    edges: &EdgeTable<T>,
    nodes: &NodeTable<T>,
) -> TablesResult<bool> {
    if edges.len() == 0 {
        return Ok(true);
    }
    let mut parent_seen = vec![0; nodes.len()];
    let mut last_parent: usize = edges[0].parent as usize;
    let mut last_child: usize = edges[0].child as usize;
    let mut last_left: i64 = edges[0].left;

    for (i, edge) in edges.iter().enumerate() {
        if edge.parent == -1 {
            // FIXME: should be name for null
            return Err(TablesError::NullParent);
        }
        if edge.child == -1 {
            // FIXME: should be name for null
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
                        if edge.left == last_left {
                            return Err(TablesError::DuplicateEdges);
                        } else if edge.left < last_left {
                            return Err(TablesError::EdgesNotSortedByLeft);
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

    return Ok(true);
}

/// A collection of node, edge, site, and mutation tables.
pub struct TableCollectionType<T> {
    length_: POSITION, // Not visible outside of this module

    pub(crate) nodes_: NodeTable<T>,
    pub(crate) edges_: EdgeTable<T>,
    pub(crate) sites_: SiteTable,
    pub(crate) mutations_: MutationTable<T>,
}

pub type TableCollection = TableCollectionType<i32>;
pub type TableCollection64 = TableCollectionType<i64>;

impl<T> TableCollectionType<T> {
    pub const fn new(genome_length: POSITION) -> TablesResult<Self> {
        if genome_length < 1 {
            return Err(TablesError::InvalidGenomeLength);
        }

        Ok(TableCollectionType {
            length_: genome_length,
            nodes_: NodeTable::<T>::new(),
            edges_: EdgeTable::<T>::new(),
            sites_: SiteTable::new(),
            mutations_: MutationTable::<T>::new(),
        })
    }

    pub fn get_length(&self) -> i64 {
        return self.length_;
    }

    /// Return immutable reference to the [mutation table](type.MutationTable.html)
    pub fn mutations(&self) -> &MutationTable<T> {
        return &self.mutations_;
    }

    /// Return immutable reference to the [edge table](type.EdgeTable.html)
    pub fn edges(&self) -> &EdgeTable<T> {
        return &self.edges_;
    }

    /// Return number of edges
    pub fn num_edges(&self) -> usize {
        return self.edges_.len();
    }

    /// Return number of nodes
    pub fn num_nodes(&self) -> usize {
        return self.nodes_.len();
    }

    /// Return immutable reference to [node table](type.NodeTable.html)
    pub fn nodes(&self) -> &NodeTable<T> {
        return &self.nodes_;
    }

    // FIXME: validate input
    pub fn node(&self, i: T) -> &Node<T> {
        return &self.nodes_[i as usize];
    }

    pub fn edge(&self, i: usize) -> &Edge<T> {
        return &self.edges_[i];
    }

    pub fn site(&self, i: usize) -> &Site {
        return &self.sites_[i];
    }

    pub fn mutation(&self, i: usize) -> &Mutation<T> {
        return &self.mutations_[i];
    }

    /// Return immutable reference to [site table](type.SiteTable.html)
    pub fn sites(&self) -> &SiteTable {
        return &self.sites_;
    }

    /// Provide an enumeration over the [node table](type.NodeTable.html)
    pub fn enumerate_nodes(&self) -> std::iter::Enumerate<std::slice::Iter<Node<T>>> {
        return self.nodes_.iter().enumerate();
    }

    /// Provide an enumeration over the [edge table](type.EdgeTable.html)
    pub fn enumerate_edges(&self) -> std::iter::Enumerate<std::slice::Iter<Edge<T>>> {
        return self.edges_.iter().enumerate();
    }

    /// Provide an enumeration over the [mutation table](type.MutationTable.html)
    pub fn enumerate_mutations(&self) -> std::iter::Enumerate<std::slice::Iter<Mutation<T>>> {
        return self.mutations_.iter().enumerate();
    }

    /// Provide an enumeration over the [site table](type.SiteTable.html)
    pub fn enumerate_sites(&self) -> std::iter::Enumerate<std::slice::Iter<Site>> {
        return self.sites_.iter().enumerate();
    }
}

pub trait TreeSequenceRecordingInterface<T> {
    // public
    fn add_node(&mut self, time: TIME, deme: T) -> TablesResult<T>;
    fn add_edge(&mut self, left: POSITION, right: POSITION, parent: T, child: T)
        -> TablesResult<T>;
    fn add_site(&mut self, position: POSITION, ancestral_state: i8) -> TablesResult<T>;
    fn add_mutation(
        &mut self,
        node: T,
        key: usize,
        site: usize,
        derived_state: i8,
        neutral: bool,
    ) -> TablesResult<T>;
    fn foo(y: Vec<i32>, x: T) -> i32 {
        return y[x];
    }
    fn sort_tables_for_simplification(&mut self) -> () {
        sort_edge_table(&self.nodes_, &mut self.edges_);
        sort_mutation_table(&self.nodes_, &mut self.edges_);
    }
}

impl TreeSequenceRecordingInterface<i32> for TableCollectionType<i32> {
    tree_sequence_recording_interface!(i32);
}

#[cfg(test)]
mod test_tables {

    use super::*;

    #[test]
    fn test_bad_genome_length() {
        let _ = TableCollection::new(0).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidGenomeLength),
            |_| assert!(false),
        );
    }

    #[test]
    fn test_add_edge() {
        let mut tables = TableCollection::new(10).unwrap();

        let result = tables.add_edge(0, 1, 2, 3).unwrap();

        assert_eq!(1, result);
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
            |_| assert!(false),
        );

        let _ = tables.add_edge(1, -1, 1, 2).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidLeftRight { found: (1, -1) }),
            |_| assert!(false),
        );
    }

    #[test]
    fn test_add_edge_bad_nodes() {
        let mut tables = TableCollection::new(10).unwrap();

        let _ = tables.add_edge(0, 1, -1, 2).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidNodeValue { found: -1 }),
            |_| assert!(false),
        );

        let _ = tables.add_edge(0, 1, 1, -2).map_or_else(
            |x: TablesError| assert_eq!(x, TablesError::InvalidNodeValue { found: -2 }),
            |_| assert!(false),
        );
    }
}
