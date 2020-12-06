use crate::tsdef::TsInt;
use thiserror::Error;

// Uses thiserror to reduce boiler plate
#[derive(Error, Debug)]
pub enum EdgeTableAddRowError {
    #[error("Bad node")]
    BadNodeValue,
    #[error("Bad position")]
    BadPosition,
}

/// A Node of a tree sequence
pub struct Node {
    /// Birth time
    pub time: i64,
    /// Population (deme) of node
    pub deme: TsInt,
}

/// An Edge is a transmission event
pub struct Edge {
    pub left: i64,
    pub right: i64,
    /// Index of parent in a [NodeTable](type.NodeTable.html)
    pub parent: TsInt,
    /// Index of child in a [NodeTable](type.NodeTable.html)
    pub child: TsInt,
}

// TODO: It would be nice to use generics here
// to allow arbitrary types for ancestral_state
// and derived_state.

/// A Site is the location and
/// ancestral state of a tables::Mutation
pub struct Site {
    pub position: i64,
    pub ancestral_state: i8,
}

/// A Mutation is the minimal information
/// needed about a mutation to track it
/// on a tree sequence.
pub struct Mutation {
    pub node: TsInt,
    pub key: usize,
    pub site: usize,
    pub derived_state: i8,
    pub neutral: bool,
}

// TODO: do these need to be pub
pub type NodeTable = Vec<Node>;
pub type EdgeTable = Vec<Edge>;
pub type SiteTable = Vec<Site>;
pub type MutationTable = Vec<Mutation>;

pub fn edge_table_add_row(
    edges: &mut EdgeTable,
    left: i64,
    right: i64,
    parent: TsInt,
    child: TsInt,
) -> Result<usize, EdgeTableAddRowError> {
    if right <= left {
        return Err(EdgeTableAddRowError::BadPosition);
    }
    if left < 0 || right < 0 {
        return Err(EdgeTableAddRowError::BadPosition);
    }
    if parent < 0 || child < 0 {
        return Err(EdgeTableAddRowError::BadNodeValue);
    }
    edges.push(Edge {
        left: left,
        right: right,
        parent: parent,
        child: child,
    });

    return Ok(edges.len());
}

// FIXME: need to validate all input params and raise errors
// if invalid.
pub fn node_table_add_row(nodes: &mut NodeTable, time: i64, deme: TsInt) -> TsInt {
    nodes.push(Node {
        time: time,
        deme: deme,
    });

    // TODO: learn if there is a way to raise error
    // automagically if overlow.
    return nodes.len() as TsInt;
}

pub fn site_table_add_row(
    sites: &mut SiteTable,
    position: i64,
    ancestral_state: i8,
) -> Result<usize, String> {
    sites.push(Site {
        position: position,
        ancestral_state: ancestral_state,
    });
    Ok(sites.len())
}

pub fn mutation_table_add_row(
    mutations: &mut MutationTable,
    node: TsInt,
    key: usize,
    site: usize,
    derived_state: i8,
    neutral: bool,
) -> usize {
    mutations.push(Mutation {
        node: node,
        key: key,
        site: site,
        derived_state: derived_state,
        neutral: neutral,
    });
    return mutations.len();
}

/// A collection of node, edge, site, and mutation tables.
pub struct TableCollection {
    length_: i64, // Not visible outside of this module

    pub(crate) nodes_: NodeTable,
    pub(crate) edges_: EdgeTable,
    pub(crate) sites_: SiteTable,
    pub(crate) mutations_: MutationTable,
}

impl TableCollection {
    pub fn add_node(&mut self, time: i64, deme: TsInt) -> TsInt {
        return node_table_add_row(&mut self.nodes_, time, deme);
    }

    /// Add an Edge
    pub fn add_edge(
        &mut self,
        left: i64,
        right: i64,
        parent: TsInt,
        child: TsInt,
    ) -> Result<usize, EdgeTableAddRowError> {
        return edge_table_add_row(&mut self.edges_, left, right, parent, child);
    }

    pub fn add_site(&mut self, position: i64, ancestral_state: i8) -> Result<usize, String> {
        if position < 0 || position >= self.length_ {
            return Err("invalid Site position".to_string());
        }
        return site_table_add_row(&mut self.sites_, position, ancestral_state);
    }

    pub fn add_mutation(
        &mut self,
        node: TsInt,
        key: usize,
        site: usize,
        derived_state: i8,
        neutral: bool,
    ) -> usize {
        return mutation_table_add_row(
            &mut self.mutations_,
            node,
            key,
            site,
            derived_state,
            neutral,
        );
    }

    pub fn get_length(&self) -> i64 {
        return self.length_;
    }

    /// Return immutable reference to the [mutation table](type.MutationTable.html)
    pub fn mutations(&self) -> &MutationTable {
        return &self.mutations_;
    }

    /// Return immutable reference to the [edge table](type.EdgeTable.html)
    pub fn edges(&self) -> &EdgeTable {
        return &self.edges_;
    }

    /// Return number of edges
    pub fn num_edges(&self) -> usize {
        return self.edges_.len();
    }

    /// Return immutable reference to [node table](type.NodeTable.html)
    pub fn nodes(&self) -> &NodeTable {
        return &self.nodes_;
    }

    /// Return immutable reference to [site table](type.SiteTable.html)
    pub fn sites(&self) -> &SiteTable {
        return &self.sites_;
    }

    // Wow, this Ord stuff takes
    // some getting used to!
    // NOTE: presumably panics if NaN/Inf show up?
    fn sort_edges(&mut self) -> () {
        // NOTE: it may by more idiomatic to
        // not use a slice here, and instead allow
        // the range-checking?
        let nslice = &self.nodes_.as_slice();
        self.edges_.sort_by(|a, b| {
            // NOTE: rust will simply NOT ALLOW
            // i32 to be an index!
            let pindex = a.parent as usize;
            let cindex = a.parent as usize;
            let ta = nslice[pindex].time;
            let tb = nslice[cindex].time;
            if ta == tb {
                if a.parent == b.parent {
                    if a.child == b.child {
                        return a.left.partial_cmp(&b.left).unwrap();
                    }
                    return a.parent.cmp(&b.parent);
                }
            }
            return ta.partial_cmp(&tb).unwrap().reverse();
        });
    }

    fn sort_mutations(&mut self) -> () {
        let sslice = &self.sites_.as_slice();
        self.mutations_.sort_by(|a, b| {
            let pa = sslice[a.site].position;
            let pb = sslice[b.site].position;
            return pa.partial_cmp(&pb).unwrap().reverse();
        });
    }

    pub fn sort_tables_for_simplification(&mut self) -> () {
        self.sort_edges();
        self.sort_mutations();
    }
}

/// Create an empty TableCollection
pub fn create_table_collection(genome_length: i64) -> TableCollection {
    return TableCollection {
        length_: genome_length,
        nodes_: NodeTable::new(),
        edges_: EdgeTable::new(),
        sites_: SiteTable::new(),
        mutations_: MutationTable::new(),
    };
}

#[cfg(test)]
mod test_tables {

    use super::*;

    #[test]
    fn test_add_edge() {
        let mut tables = create_table_collection(10);

        let result = tables.add_edge(0, 1, 2, 3);

        match result {
            Ok(_usize) => assert!(true),
            Err(EdgeTableAddRowError::BadPosition) => assert!(false),
            Err(EdgeTableAddRowError::BadNodeValue) => assert!(false),
        }

        assert_eq!(1, tables.edges().len());
        assert_eq!(1, tables.num_edges());
    }

    #[test]
    fn test_add_bad_edges() {
        let mut tables = create_table_collection(10);

        // Add bad coordinates
        let result = tables.add_edge(1, -1, 1, 2);
        match result {
            Ok(_usize) => assert!(false),
            Err(EdgeTableAddRowError::BadPosition) => assert!(true),
            Err(EdgeTableAddRowError::BadNodeValue) => assert!(false),
        }

        let result = tables.add_edge(-1, 2, 1, 2);
        match result {
            Ok(_usize) => assert!(false),
            Err(EdgeTableAddRowError::BadPosition) => assert!(true),
            Err(EdgeTableAddRowError::BadNodeValue) => assert!(false),
        }

        // Negative child/parent values
        let result = tables.add_edge(0, 1, 1, -2);
        match result {
            Ok(_usize) => assert!(false),
            Err(EdgeTableAddRowError::BadPosition) => assert!(false),
            Err(EdgeTableAddRowError::BadNodeValue) => assert!(true),
        }

        let result = tables.add_edge(0, 1, -1, 2);
        match result {
            Ok(_usize) => assert!(false),
            Err(EdgeTableAddRowError::BadPosition) => assert!(false),
            Err(EdgeTableAddRowError::BadNodeValue) => assert!(true),
        }
    }
}
