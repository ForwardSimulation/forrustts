use crate::tsdef::{IdType, Position, Time, NULL_ID};
use std::cmp::Ordering;
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum TablesError {
    #[error("Invalid genome length")]
    InvalidGenomeLength,
    #[error("Invalid node: {found:?}")]
    InvalidNodeValue { found: IdType },
    #[error("Invalid value for position: {found:?}")]
    InvalidPosition { found: Position },
    #[error("Invalid position range: {found:?}")]
    InvalidLeftRight { found: (Position, Position) },
    #[error("Invalid value for time: {found:?}")]
    InvalidTime { found: Time },
    #[error("Invalid value for deme: {found:?}")]
    InvalidDeme { found: IdType },
    #[error("Parent is NULL_ID")]
    NullParent,
    #[error("Child is NULL_ID")]
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
pub struct Node {
    /// Birth time
    pub time: Time,
    /// Population (deme) of node
    pub deme: IdType,
}

/// An Edge is a transmission event
#[derive(Copy, Clone)]
pub struct Edge {
    pub left: Position,
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
pub struct Site {
    pub position: Position,
    pub ancestral_state: Option<Vec<u8>>,
}

/// A MutationRecord is the minimal information
/// needed about a mutation to track it
/// on a tree sequence.
pub struct MutationRecord {
    pub node: IdType,
    pub key: usize,
    pub site: usize,
    pub derived_state: Option<Vec<u8>>,
    pub neutral: bool,
}

// TODO: do these need to be pub
pub type NodeTable = Vec<Node>;
pub type EdgeTable = Vec<Edge>;
pub type SiteTable = Vec<Site>;
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

fn time_non_negative(x: Time) -> TablesResult<()> {
    if x < 0 {
        Err(TablesError::InvalidTime { found: x })
    } else {
        Ok(())
    }
}

fn deme_non_negative(x: IdType) -> TablesResult<()> {
    if x < 0 {
        Err(TablesError::InvalidDeme { found: x })
    } else {
        Ok(())
    }
}

pub fn edge_table_add_row(
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

    Ok(edges.len() as IdType)
}

// FIXME: need to validate all input params and raise errors
// if invalid.
pub fn node_table_add_row(nodes: &mut NodeTable, time: Time, deme: IdType) -> TablesResult<IdType> {
    time_non_negative(time)?;
    deme_non_negative(deme)?;
    nodes.push(Node { time, deme });

    Ok((nodes.len() - 1) as IdType)
}

pub fn site_table_add_row(
    sites: &mut SiteTable,
    position: Position,
    ancestral_state: Option<Vec<u8>>,
) -> TablesResult<IdType> {
    position_non_negative(position)?;
    sites.push(Site {
        position,
        ancestral_state,
    });
    Ok(sites.len() as IdType)
}

pub fn mutation_table_add_row(
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
    Ok(mutations.len() as IdType)
}

fn sort_edge_table(nodes: &[Node], edges: &mut EdgeTable) {
    // NOTE: it may by more idiomatic to
    // not use a slice here, and instead allow
    // the range-checking?
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

fn sort_mutation_table(sites: &[Site], mutations: &mut MutationTable) {
    mutations.sort_by(|a, b| {
        let pa = sites[a.site].position;
        let pb = sites[b.site].position;
        pa.cmp(&pb).reverse()
    });
}

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
pub struct TableCollection {
    length_: Position, // Not visible outside of this module

    pub(crate) nodes_: NodeTable,
    pub(crate) edges_: EdgeTable,
    pub(crate) sites_: SiteTable,
    pub(crate) mutations_: MutationTable,
}

impl TableCollection {
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

    pub fn add_node(&mut self, time: Time, deme: IdType) -> TablesResult<IdType> {
        node_table_add_row(&mut self.nodes_, time, deme)
    }

    /// Add an Edge
    pub fn add_edge(
        &mut self,
        left: Position,
        right: Position,
        parent: IdType,
        child: IdType,
    ) -> TablesResult<IdType> {
        edge_table_add_row(&mut self.edges_, left, right, parent, child)
    }

    pub fn add_site(
        &mut self,
        position: Position,
        ancestral_state: Option<Vec<u8>>,
    ) -> TablesResult<IdType> {
        if position >= self.length_ {
            return Err(TablesError::InvalidPosition { found: position });
        }
        site_table_add_row(&mut self.sites_, position, ancestral_state)
    }

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

    pub fn get_length(&self) -> Position {
        self.length_
    }

    /// Return immutable reference to the [mutation table](type.MutationTable.html)
    pub fn mutations(&self) -> &MutationTable {
        &self.mutations_
    }

    /// Return immutable reference to the [edge table](type.EdgeTable.html)
    pub fn edges(&self) -> &EdgeTable {
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
    pub fn nodes(&self) -> &NodeTable {
        &self.nodes_
    }

    pub fn node(&self, i: IdType) -> &Node {
        &self.nodes_[i as usize]
    }

    pub fn edge(&self, i: IdType) -> &Edge {
        &self.edges_[i as usize]
    }

    pub fn site(&self, i: IdType) -> &Site {
        &self.sites_[i as usize]
    }

    pub fn mutation(&self, i: IdType) -> &MutationRecord {
        &self.mutations_[i as usize]
    }

    /// Return immutable reference to [site table](type.SiteTable.html)
    pub fn sites(&self) -> &SiteTable {
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

    pub fn sort_tables_for_simplification(&mut self) {
        sort_edge_table(&self.nodes_, &mut self.edges_);
        sort_mutation_table(&self.sites_, &mut self.mutations_);
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
}
