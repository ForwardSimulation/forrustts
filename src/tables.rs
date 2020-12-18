use crate::tsdef::{POSITION, TIME};
use paste;
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
    InvalidDeme { found: i64 },
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

fn time_non_negative(x: TIME) -> TablesResult<()> {
    if x < 0 {
        return Err(TablesError::InvalidTime { found: x });
    }
    Ok(())
}

macro_rules! node_non_negative {
    ($itype: ty) => {
        fn node_non_negative(x: $itype) -> TablesResult<()> {
            if x < 0 {
                return Err(TablesError::InvalidNodeValue { found: x as i64 });
            }
            Ok(())
        }
    };
}

macro_rules! deme_non_negative {
    ($itype: ty) => {
        fn deme_non_negative(x: $itype) -> TablesResult<()> {
            if x < 0 {
                return Err(TablesError::InvalidDeme { found: x as i64 });
            }
            Ok(())
        }
    };
}

macro_rules! add_edge {
    ($itype: ty) => {
        fn add_edge(
            //edges: &mut EdgeTable<$itype>,
            &mut self,
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
            Self::node_non_negative(parent)?;
            Self::node_non_negative(child)?;

            self.edges_.push(Edge::<$itype> {
                left: left,
                right: right,
                parent: parent,
                child: child,
            });

            Ok(self.edges_.len() as $itype)
        }
    };
}

macro_rules! add_node {
    ($itype: ty) => {
        fn add_node(
            //nodes: &mut NodeTable<$itype>,
            &mut self,
            time: TIME,
            deme: $itype,
        ) -> TablesResult<$itype> {
            time_non_negative(time)?;
            Self::deme_non_negative(deme)?;
            self.nodes_.push(Node::<$itype> {
                time: time,
                deme: deme,
            });

            // TODO: learn if there is a way to raise error
            // automagically if overlow.
            Ok((self.nodes_.len() - 1) as $itype)
        }
    };
}

macro_rules! add_site {
    ($itype: ty) => {
        fn add_site(
            //sites: &mut SiteTable,
            &mut self,
            position: POSITION,
            ancestral_state: i8,
        ) -> TablesResult<$itype> {
            position_non_negative(position)?;
            self.sites_.push(Site {
                position: position,
                ancestral_state: ancestral_state,
            });
            Ok(self.sites_.len() as $itype)
        }
    };
}

macro_rules! add_mutation {
    ($itype: ty) => {
        fn add_mutation(
            //mutations: &mut MutationTable<$itype>,
            &mut self,
            node: $itype,
            key: usize,
            site: usize,
            derived_state: i8,
            neutral: bool,
        ) -> TablesResult<$itype> {
            Self::node_non_negative(node)?;
            self.mutations_.push(Mutation::<$itype> {
                node: node,
                key: key,
                site: site,
                derived_state: derived_state,
                neutral: neutral,
            });
            Ok(self.mutations_.len() as $itype)
        }
    };
}

macro_rules! sort_edge_table {
    ($itype: ty, $fn: expr) => {
        paste::item! {
        fn $fn(nodes: &NodeTable<$itype>, edges: &mut EdgeTable<$itype>) -> () {
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
        }
    };
}

macro_rules! sort_mutation_table {
    ($itype: ty, $fn: expr) => {
        paste::item! {
        fn $fn(sites: &SiteTable, mutations: &mut MutationTable<$itype>) -> () {
            let sslice = &sites.as_slice();
            mutations.sort_by(|a, b| {
                let pa = sslice[a.site].position;
                let pb = sslice[b.site].position;
                return pa.partial_cmp(&pb).unwrap().reverse();
            });
        }
        }
    };
}

macro_rules! sort_tables_for_simplification {
    ($itype: ty, $edge_fn: expr, $mut_fn: expr) => {
        fn sort_tables_for_simplification(&mut self) -> () {
            $edge_fn(&self.nodes_, &mut self.edges_);
            $mut_fn(&self.sites_, &mut self.mutations_);
        }
    };
}

macro_rules! auxilliary_sorting_functions {
    ($itype:ty, $edge_fn: expr, $mut_fn: expr) => {
        sort_mutation_table!($itype, $mut_fn);
        sort_edge_table!($itype, $edge_fn);
    };
}

macro_rules! tree_sequence_recording_interface {
    ($itype: ty, $edge_fn: expr, $mut_fn: expr) => {
        // NOTE: should be hidden
        node_non_negative!($itype);
        deme_non_negative!($itype);
        // The public interfact
        add_node!($itype);
        add_edge!($itype);
        add_site!($itype);
        add_mutation!($itype);
        sort_tables_for_simplification!($itype, $edge_fn, $mut_fn);
        validate_edge_table!($itype);
    };
}

macro_rules! validate_edge_table {
    ($itype: ty) => {
        fn validate_edge_table(&self) -> TablesResult<bool> {
            let edges = self.edges_.as_slice();
            let nodes = self.nodes_.as_slice();
            let len = self.get_length();
            if edges.len() == 0 {
                return Ok(true);
            }
            let mut parent_seen = vec![0; nodes.len()];
            let mut last_parent: usize = edges[0].parent as usize;
            let mut last_child: usize = edges[0].child as usize;
            let mut last_left: i64 = edges[0].left;

            for (i, edge) in edges.iter().enumerate() {
                if edge.parent == Self::NULL_ID {
                    return Err(TablesError::NullParent);
                }
                if edge.child == Self::NULL_ID {
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
    };
}

/// A collection of node, edge, site, and mutation tables.
pub struct TableCollectionType<T> {
    length_: POSITION, // Not visible outside of this module

    pub(crate) nodes_: NodeTable<T>,
    pub(crate) edges_: EdgeTable<T>,
    pub(crate) sites_: SiteTable,
    pub(crate) mutations_: MutationTable<T>,
}

/// The current "canonical" tables with 32 bit identifiers
pub type TableCollection = TableCollectionType<i32>;
/// Tables with 64 bit identifiers
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

    // TODO: all of the below should be taking in ID types:

    // FIXME: validate input
    pub fn node(&self, i: usize) -> &Node<T> {
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
    type IdType;
    const NULL_ID: Self::IdType;
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
    fn sort_tables_for_simplification(&mut self) -> ();
    fn validate_edge_table(&self) -> TablesResult<bool>;

    // Trait interfaces are all "public",
    // so we hide some details from the docs.
    #[doc(hidden)]
    fn node_non_negative(x: T) -> TablesResult<()>;
    #[doc(hidden)]
    fn deme_non_negative(x: T) -> TablesResult<()>;
}

// TODO: explore simplifying these macros by making
// the internal sorting functions trait members
// w/hidden docs?

auxilliary_sorting_functions!(i32, sort_edge_table_i32, sort_mutation_table_i32);

impl TreeSequenceRecordingInterface<i32> for TableCollectionType<i32> {
    type IdType = i32;
    const NULL_ID: Self::IdType = -1;
    tree_sequence_recording_interface!(i32, sort_edge_table_i32, sort_mutation_table_i32);
}

auxilliary_sorting_functions!(i64, sort_edge_table_i64, sort_mutation_table_i64);

impl TreeSequenceRecordingInterface<i64> for TableCollectionType<i64> {
    type IdType = i64;
    const NULL_ID: Self::IdType = -1;
    tree_sequence_recording_interface!(i64, sort_edge_table_i64, sort_mutation_table_i64);
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
