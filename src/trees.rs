use crate::{IdType, IdVec, Position, NULL_ID};
use bitflags::bitflags;

bitflags! {
    pub struct TreeFlags: u32 {
        const NONE = 0;
        const TRACK_SAMPLES = 1 << 0;
    }
}

#[derive(Copy, Clone)]
pub struct TopologyData {
    parent: IdType,
    left_child: IdType,
    right_child: IdType,
    left_sib: IdType,
    right_sib: IdType,
    left_sample: IdType,
    right_sample: IdType,
    next_sample: IdType,
    leaf_counts: IdType,
}

impl Default for TopologyData {
    fn default() -> Self {
        Self {
            parent: NULL_ID,
            left_child: NULL_ID,
            right_child: NULL_ID,
            left_sib: NULL_ID,
            right_sib: NULL_ID,
            left_sample: NULL_ID,
            right_sample: NULL_ID,
            next_sample: NULL_ID,
            leaf_counts: 0,
        }
    }
}

pub struct Tree<'treeseq> {
    topology: Vec<TopologyData>,
    left_root: IdType,
    above_sample: Vec<i8>,
    left: Position,
    right: Position,
    samples: IdVec,
    sample_index_map: IdVec, // TODO: decide if this is better as usize.
    flags: TreeFlags,
    treeseq: &'treeseq TreeSequence<'treeseq>,
    // The following help implement StreamingIterator
    input_edge_index: usize,
    output_edge_index: usize,
    x: Position,
    advanced: bool,
}

impl<'treeseq> Tree<'treeseq> {
    fn new_internal(treeseq: &'treeseq TreeSequence, samples: &[IdType], flags: TreeFlags) -> Self {
        Self {
            topology: vec![TopologyData::default(); treeseq.tables.num_nodes()],
            left_root: NULL_ID,
            above_sample: vec![0; treeseq.tables.num_nodes()],
            left: Position::MIN,
            right: Position::MIN,
            samples: samples.to_vec(),
            sample_index_map: vec![NULL_ID; treeseq.tables.num_nodes()],
            flags,
            treeseq,
            input_edge_index: 0,
            output_edge_index: 0,
            x: 0,
            advanced: false,
        }
    }

    fn init_samples(&mut self) {
        for (i, s) in self.samples.iter().enumerate() {
            if self.sample_index_map[*s as usize] != NULL_ID {
                panic!("invalid sample list");
            }
            self.sample_index_map[*s as usize] = i as IdType;
            if let Some(row) = self.topology.get_mut(*s as usize) {
                row.left_sample = self.sample_index_map[*s as usize];
                row.right_sample = self.sample_index_map[*s as usize];
                row.leaf_counts = 1;
                self.above_sample[*s as usize] = 1;

                // Initialize roots
                if i < self.samples.len() - 1 {
                    row.right_sib = *unsafe { self.samples.get_unchecked(i + 1) };
                }
                if i > 0 {
                    row.left_sib = *unsafe { self.samples.get_unchecked(i - 1) };
                }
            } else {
                panic!("expected Some(mut row)");
            }
        }
    }

    fn update_incoming_leaf_count(&mut self, parent: IdType, child: IdType) {
        let mut u = parent;
        let lc = self.topology[child as usize].leaf_counts;
        if lc == 0 {
            return;
        }
        while u != NULL_ID {
            self.topology[u as usize].leaf_counts += lc;
            u = self.topology[u as usize].parent;
        }
    }

    fn update_incoming_roots(&mut self, parent: IdType, child: IdType, lsib: IdType, rsib: IdType) {
        if self.above_sample[child as usize] > 0 {
            let mut x = parent;
            let mut root = x;
            let mut above_sample = false;

            while x != NULL_ID && !above_sample {
                above_sample = self.above_sample[x as usize] > 0;
                // c is above_sample and p is c's parent.
                // Thus, all parents to p are above_sample, too.
                self.above_sample[x as usize] = 1;
                root = x;
                x = self.topology[x as usize].parent;
            }

            if !above_sample {
                // If we are here, then the above loop terminated
                // by encountering a NULL node, because above_sample[x]
                // must have been 0 for all x. However, because c is
                // above sample, all nodes encountered have been update
                // to be above_sample as well. Thus, the new value of root
                // replaces c in the root list.

                if lsib != NULL_ID {
                    self.topology[lsib as usize].right_sib = root;
                }
                if rsib != NULL_ID {
                    self.topology[rsib as usize].left_sib = root;
                }
                self.topology[root as usize].left_sib = lsib;
                self.topology[root as usize].right_sib = lsib;
                self.left_root = root;
            } else {
                // If we are here, then we encountered a node
                // ancestral to c where above_sample == 1.
                // Thus, c can no longer be a root.  If the current
                // p is also a c in a later call to this function, then
                // it may also be removed, etc..
                self.left_root = NULL_ID;
                if lsib != NULL_ID {
                    self.topology[lsib as usize].right_sib = rsib;
                    self.left_root = lsib;
                }
                if rsib != NULL_ID {
                    self.topology[rsib as usize].left_sib = lsib;
                    self.left_root = rsib;
                }
            }
        }
    }

    fn update_outgoing_leaf_count(&mut self, parent: IdType, child: IdType) {
        let mut u = parent;
        let lc = self.topology[child as usize].leaf_counts;
        if lc == 0 {
            return;
        }
        while u != NULL_ID {
            self.topology[u as usize].leaf_counts -= lc;
            u = self.topology[u as usize].parent;
        }
    }

    fn update_outgoing_roots(&mut self, parent: IdType, child: IdType) {
        if self.above_sample[child as usize] == 1 {
            let mut x = parent;
            let mut root = x;
            let mut above_sample = false;

            while x != NULL_ID && !above_sample {
                above_sample = self.sample_index_map[x as usize] != NULL_ID;
                let mut lc = self.topology[x as usize].left_child;
                while lc != NULL_ID && !above_sample {
                    above_sample = above_sample || self.above_sample[lc as usize] > 0;
                    lc = self.topology[lc as usize].left_sib;
                }
                if above_sample {
                    self.above_sample[x as usize] = 1;
                }
                root = x;
                x = self.topology[x as usize].parent;
            }

            // Now, root refers to the most ancient
            // ancestor of parent found in the above loop
            if !above_sample {
                // remove root from list of roots
                let lroot = self.topology[root as usize].left_sib;
                let rroot = self.topology[root as usize].right_sib;
                self.left_root = NULL_ID;
                if lroot != NULL_ID {
                    self.topology[lroot as usize].right_sib = rroot;
                    self.left_root = lroot;
                }
                if rroot != NULL_ID {
                    self.topology[rroot as usize].left_sib = lroot;
                    self.left_root = rroot;
                }
                self.topology[root as usize].left_sib = NULL_ID;
                self.topology[root as usize].right_sib = NULL_ID;
            }
            if self.left_root != NULL_ID {
                let lroot = self.topology[self.left_root as usize].left_sib;
                if lroot != NULL_ID {
                    self.topology[lroot as usize].right_sib = child;
                }
                self.topology[child as usize].left_sib = lroot;
                self.topology[self.left_root as usize].left_sib = child;
            }
            self.topology[child as usize].right_sib = self.left_root;
            self.left_root = child;
        }
    }

    fn update_samples_list(&mut self, node: IdType) {
        assert!(self.flags.contains(TreeFlags::TRACK_SAMPLES));

        let sample_map = self.sample_index_map.as_slice();
        let topo = self.topology.as_mut_slice();
        let mut n = node;

        while n != NULL_ID {
            let sample_index = sample_map[n as usize];
            if sample_index != NULL_ID {
                topo[n as usize].right_sample = topo[n as usize].left_sample;
            } else {
                topo[n as usize].left_sample = NULL_ID;
                topo[n as usize].right_sample = NULL_ID;
            }

            let mut v = topo[n as usize].left_child;
            while v != NULL_ID {
                if topo[v as usize].left_sample != NULL_ID {
                    assert!(topo[v as usize].right_sample != NULL_ID);
                    if topo[n as usize].left_sample == NULL_ID {
                        topo[n as usize].left_sample = topo[v as usize].left_sample;
                        topo[n as usize].right_sample = topo[v as usize].right_sample;
                    } else {
                        let nright = topo[n as usize].right_sample as usize;
                        let vleft = topo[v as usize].left_sample;
                        topo[nright].next_sample = vleft;
                        topo[n as usize].right_sample = topo[v as usize].right_sample;
                    }
                }
                v = topo[v as usize].right_sib;
            }
            n = topo[n as usize].parent;
        }
    }

    pub fn new(treeseq: &'treeseq TreeSequence, samples: &[IdType], flags: TreeFlags) -> Self {
        assert!(!samples.is_empty());
        let mut rv = Self::new_internal(treeseq, samples, flags);
        rv.init_samples();
        rv.left_root = rv.samples[0];
        rv
    }
}

/// Left-to-right iteration of trees.
impl<'treeseq> streaming_iterator::StreamingIterator for Tree<'treeseq> {
    type Item = Tree<'treeseq>;

    fn advance(&mut self) {
        let tables = &self.treeseq.tables;
        let edge_table = self.treeseq.tables.edges_.as_slice();
        let edge_input_order = tables.edge_input_order.as_slice();
        let edge_output_order = tables.edge_output_order.as_slice();
        if self.input_edge_index < edge_input_order.len() || self.x < tables.genome_length() {
            while self.output_edge_index < edge_output_order.len()
                && edge_table[edge_output_order[self.output_edge_index]].right == self.x
            {
                let current_edge = edge_table[edge_output_order[self.output_edge_index]];
                let lsib = self.topology[current_edge.child as usize].left_sib;
                let rsib = self.topology[current_edge.child as usize].right_sib;

                if lsib == NULL_ID {
                    self.topology[current_edge.parent as usize].left_child = rsib;
                } else {
                    self.topology[lsib as usize].right_sib = rsib;
                }
                if rsib == NULL_ID {
                    self.topology[current_edge.parent as usize].right_child = lsib;
                } else {
                    self.topology[rsib as usize].left_sib = lsib;
                }
                let child_topo = &mut self.topology[current_edge.child as usize];
                child_topo.parent = NULL_ID;
                child_topo.left_sib = NULL_ID;
                child_topo.right_sib = NULL_ID;

                self.update_outgoing_leaf_count(current_edge.parent, current_edge.child);

                if self.flags.contains(TreeFlags::TRACK_SAMPLES) {
                    self.update_samples_list(current_edge.parent);
                }
                self.update_outgoing_roots(current_edge.parent, current_edge.child);
                self.output_edge_index += 1;
            }
            while self.input_edge_index < edge_input_order.len()
                && edge_table[edge_input_order[self.input_edge_index]].left == self.x
            {
                let current_edge = edge_table[edge_input_order[self.input_edge_index]];
                let rchild = self.topology[current_edge.parent as usize].right_child;
                let lsib = self.topology[current_edge.child as usize].left_sib;
                let rsib = self.topology[current_edge.child as usize].right_sib;

                if rchild == NULL_ID {
                    self.topology[current_edge.parent as usize].left_child = current_edge.child;
                    self.topology[current_edge.child as usize].left_sib = NULL_ID;
                    self.topology[current_edge.child as usize].right_sib = NULL_ID;
                } else {
                    self.topology[rchild as usize].right_sib = current_edge.child;
                    self.topology[current_edge.child as usize].left_sib = rchild;
                    self.topology[current_edge.child as usize].right_sib = NULL_ID;
                }
                self.topology[current_edge.child as usize].parent = current_edge.parent;
                self.topology[current_edge.parent as usize].right_child = current_edge.child;

                self.update_incoming_leaf_count(current_edge.parent, current_edge.child);
                if self.flags.contains(TreeFlags::TRACK_SAMPLES) {
                    self.update_samples_list(current_edge.parent);
                }
                self.update_incoming_roots(current_edge.parent, current_edge.child, lsib, rsib);
                self.input_edge_index += 1;
            }

            // This is a big "gotcha".
            // The root tracking functions will sometimes
            // result in left_root not actually being the left_root.
            // We loop through the left_sibs to fix that.
            if self.left_root != NULL_ID {
                while self.topology[self.left_root as usize].left_sib != NULL_ID {
                    self.left_root = self.topology[self.left_root as usize].left_sib;
                }
            }

            let mut right = tables.genome_length();
            if self.input_edge_index < edge_input_order.len() {
                right = std::cmp::min(
                    right,
                    edge_table[edge_input_order[self.input_edge_index]].left,
                );
            }
            if self.output_edge_index < edge_output_order.len() {
                right = std::cmp::min(
                    right,
                    edge_table[edge_output_order[self.output_edge_index]].right,
                );
            }
            self.left = self.x;
            self.right = right;
            self.x = right;
            self.advanced = true;
        } else {
            self.advanced = false;
        }
    }

    fn get(&self) -> Option<&Self::Item> {
        match self.advanced {
            true => Some(&self),
            false => None,
        }
    }
}

/// Error type related to [``TreeSequence``] and [``Tree``].
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum TreesError {
    /// Raised by [``TreeSequence::new``].
    #[error("Tables not indexed.")]
    TablesNotIndexed,
}

pub struct TreeSequence<'a> {
    tables: &'a crate::TableCollection,
}

/// Result type for operations on trees and tree sequences.
pub type TreesResult<T> = Result<T, TreesError>;

impl<'a> TreeSequence<'a> {
    pub fn new(tables: &'a crate::TableCollection) -> TreesResult<Self> {
        Ok(Self { tables })
    }

    pub fn tables(&self) -> &'a crate::TableCollection {
        self.tables
    }

    pub fn tree_iter(&self, samples: &[IdType], flags: TreeFlags) -> Tree<'_> {
        Tree::new(self, samples, flags)
    }
}

#[cfg(test)]
mod test_trees {

    #[test]
    fn test_treeseq_creation_and_table_access() {
        let mut tables = crate::TableCollection::new(100).unwrap();
        tables.add_edge(0, 1, 0, 1).unwrap();

        let ts = super::TreeSequence::new(&tables).unwrap();

        let tref = ts.tables();
        assert_eq!(tref.edges().len(), 1);
    }

    #[test]
    fn test_treeseq_creation_and_tree_creation() {
        let mut tables = crate::TableCollection::new(100).unwrap();
        tables.add_edge(0, 1, 0, 1).unwrap();
        tables.add_node(0, 0).unwrap();
        tables.add_node(1, 0).unwrap();

        let ts = super::TreeSequence::new(&tables).unwrap();
        let _ = ts.tree_iter(&[0], super::TreeFlags::empty());
    }

    fn make_small_table_collection_two_trees() -> crate::TableCollection {
        // The two trees are:
        //  0
        // +++
        // | |  1
        // | | +++
        // 2 3 4 5

        //     0
        //   +-+-+
        //   1   |
        // +-+-+ |
        // 2 4 5 3

        let mut tables = crate::TableCollection::new(1000).unwrap();
        tables.add_node(0, 0).unwrap();
        tables.add_node(1, 0).unwrap();
        tables.add_node(2, 0).unwrap(); // nodes 2-5 are samples
        tables.add_node(2, 0).unwrap();
        tables.add_node(2, 0).unwrap();
        tables.add_node(2, 0).unwrap();
        tables.add_edge(500, 1000, 0, 1).unwrap();
        tables.add_edge(0, 500, 0, 2).unwrap();
        tables.add_edge(0, 1000, 0, 3).unwrap();
        tables.add_edge(500, 1000, 1, 2).unwrap();
        tables.add_edge(0, 1000, 1, 4).unwrap();
        tables.add_edge(0, 1000, 1, 5).unwrap();
        tables.sort_tables(crate::TableSortingFlags::default());
        tables
            .validate(crate::TableValidationFlags::VALIDATE_ALL)
            .unwrap();
        tables
            .build_indexes(crate::IndexTablesFlags::empty())
            .unwrap();
        tables
    }

    #[test]
    fn test_two_trees() {
        use streaming_iterator::StreamingIterator;
        let tables = make_small_table_collection_two_trees();
        let treeseq = super::TreeSequence::new(&tables).unwrap();
        let samples: crate::IdVec = vec![2, 3, 4, 5];

        for i in &samples {
            assert_eq!(tables.node(*i).time, 2);
        }

        let mut tree_iter = treeseq.tree_iter(&samples, super::TreeFlags::TRACK_SAMPLES);
        let mut ntrees = 0;
        while tree_iter.next().is_some() {
            ntrees += 1;
        }
        assert_eq!(ntrees, 2);
    }
}
