use crate::newtypes::{NodeId, Position, Time};
use bitflags::bitflags;

bitflags! {
    /// Modify the behavior of [`TreeSequence::tree_iterator`].
    #[derive(Default)]
    pub struct TreeFlags: u32 {
        /// Keep track of which sample nodes descend from
        /// each node.  This tracking is relatively expensive.
        const TRACK_SAMPLES = 1 << 0;
    }
}

/// Data describing the toplological relationship
/// between [`NodeId`] in a [`Tree`].
///
/// For a [`TreeSequence`] whose tables have `n`
/// nodes, there are `n` instances of this
/// struct.
///
/// For a given instance, the fields provide
/// the id of other nodes of specific relationships
/// in the same tree.
///
/// Some fields may be equal to [`NodeId::NULL`],
/// indicating that the current instance is a root
/// or leaf node, for example.
#[derive(Copy, Clone)]
struct TopologyData {
    parent: NodeId,
    left_child: NodeId,
    right_child: NodeId,
    left_sib: NodeId,
    right_sib: NodeId,
    left_sample: NodeId,
    right_sample: NodeId,
    next_sample: NodeId,
    leaf_counts: i32,
}

impl Default for TopologyData {
    fn default() -> Self {
        Self {
            parent: NodeId::NULL,
            left_child: NodeId::NULL,
            right_child: NodeId::NULL,
            left_sib: NodeId::NULL,
            right_sib: NodeId::NULL,
            left_sample: NodeId::NULL,
            right_sample: NodeId::NULL,
            next_sample: NodeId::NULL,
            leaf_counts: 0,
        }
    }
}

trait NodeIterator {
    fn next_node(&mut self);
    fn current_node(&mut self) -> Option<NodeId>;
}

/// Specify the traversal order used by
/// [`Tree::traverse_nodes`].
pub enum NodeTraversalOrder {
    ///Preorder traversal, starting at the root(s) of a [`Tree`].
    ///For trees with multiple roots, start at the left root,
    ///traverse to tips, proceeed to the next root, etc..
    Preorder,
}

struct PreorderNodeIterator<'a> {
    node_stack: Vec<NodeId>,
    tree: &'a Tree<'a>,
    current_node_: Option<NodeId>,
}

impl<'a> PreorderNodeIterator<'a> {
    fn new(tree: &'a Tree) -> Self {
        let mut rv = PreorderNodeIterator {
            node_stack: tree.roots_to_vec(),
            tree,
            current_node_: None,
        };
        rv.node_stack.reverse();
        rv
    }
}

impl NodeIterator for PreorderNodeIterator<'_> {
    fn next_node(&mut self) {
        self.current_node_ = self.node_stack.pop();
        if let Some(u) = self.current_node_ {
            let mut c = self.tree.right_child(u).unwrap();
            while c != NodeId::NULL {
                self.node_stack.push(c);
                c = self.tree.left_sib(c).unwrap();
            }
        };
    }

    fn current_node(&mut self) -> Option<NodeId> {
        self.current_node_
    }
}

iterator_for_nodeiterator!(PreorderNodeIterator<'_>);

struct RootIterator<'a> {
    current_root: Option<NodeId>,
    next_root: NodeId,
    tree: &'a Tree<'a>,
}

impl<'a> RootIterator<'a> {
    fn new(tree: &'a Tree) -> Self {
        RootIterator {
            current_root: None,
            next_root: tree.left_root,
            tree,
        }
    }
}

impl NodeIterator for RootIterator<'_> {
    fn next_node(&mut self) {
        self.current_root = match self.next_root {
            NodeId::NULL => None,
            r => {
                assert!(r >= 0);
                let cr = Some(r);
                self.next_root = self.tree.right_sib(r).unwrap();
                cr
            }
        };
    }

    fn current_node(&mut self) -> Option<NodeId> {
        self.current_root
    }
}

iterator_for_nodeiterator!(RootIterator<'_>);

struct ChildIterator<'a> {
    current_child: Option<NodeId>,
    next_child: NodeId,
    tree: &'a Tree<'a>,
}

impl<'a> ChildIterator<'a> {
    fn new(tree: &'a Tree, u: NodeId) -> Self {
        let c = tree.left_child(u).unwrap();

        ChildIterator {
            current_child: None,
            next_child: c,
            tree,
        }
    }
}

impl NodeIterator for ChildIterator<'_> {
    fn next_node(&mut self) {
        self.current_child = match self.next_child {
            NodeId::NULL => None,
            r => {
                assert!(r >= 0);
                let cr = Some(r);
                self.next_child = self.tree.right_sib(r).unwrap();
                cr
            }
        };
    }

    fn current_node(&mut self) -> Option<NodeId> {
        self.current_child
    }
}

iterator_for_nodeiterator!(ChildIterator<'_>);

struct ParentsIterator<'a> {
    current_node: Option<NodeId>,
    next_node: NodeId,
    tree: &'a Tree<'a>,
}

impl<'a> ParentsIterator<'a> {
    fn new(tree: &'a Tree, u: NodeId) -> Self {
        ParentsIterator {
            current_node: None,
            next_node: u,
            tree,
        }
    }
}

impl NodeIterator for ParentsIterator<'_> {
    fn next_node(&mut self) {
        self.current_node = match self.next_node {
            NodeId::NULL => None,
            r => {
                assert!(r >= 0);
                let cr = Some(r);
                self.next_node = self.tree.parent(r).unwrap();
                cr
            }
        };
    }

    fn current_node(&mut self) -> Option<NodeId> {
        self.current_node
    }
}

iterator_for_nodeiterator!(ParentsIterator<'_>);

struct SamplesIterator<'a> {
    current_node: Option<NodeId>,
    next_sample_index: NodeId,
    last_sample_index: NodeId,
    tree: &'a Tree<'a>,
}

impl<'a> SamplesIterator<'a> {
    fn new(tree: &'a Tree, u: NodeId) -> Self {
        SamplesIterator {
            current_node: None,
            next_sample_index: tree.left_sample(u).unwrap(),
            last_sample_index: tree.right_sample(u).unwrap(),
            tree,
        }
    }
}

impl NodeIterator for SamplesIterator<'_> {
    fn next_node(&mut self) {
        self.current_node = match self.next_sample_index {
            NodeId::NULL => None,
            r => {
                if r == self.last_sample_index {
                    let cr = Some(self.tree.samples[r.0 as usize]);
                    self.next_sample_index = NodeId::NULL;
                    cr
                } else {
                    assert!(r >= 0);
                    let cr = Some(self.tree.samples[r.0 as usize]);
                    self.next_sample_index = self.tree.topology[r.0 as usize].next_sample;
                    cr
                }
            }
        };
    }

    fn current_node(&mut self) -> Option<NodeId> {
        self.current_node
    }
}

iterator_for_nodeiterator!(SamplesIterator<'_>);

/// A tree is the genealogy of a non-recombining
/// segment of a genome.  A [`TreeSequence`] contains
/// the information needed to efficiently build trees
/// and iterate over each tree in a genome.
pub struct Tree<'treeseq> {
    topology: Vec<TopologyData>,
    left_root: NodeId,
    above_sample: Vec<i8>,
    left: Position,
    right: Position,
    samples: &'treeseq [NodeId],
    sample_index_map: Vec<NodeId>, // TODO: decide if this is better as usize.
    flags: TreeFlags,
    treeseq: &'treeseq TreeSequence,
    // The following help implement StreamingIterator
    input_edge_index: usize,
    output_edge_index: usize,
    x: Position,
    advanced: bool,
}

impl<'treeseq> Tree<'treeseq> {
    fn new_internal(treeseq: &'treeseq TreeSequence, flags: TreeFlags) -> Self {
        Self {
            topology: vec![TopologyData::default(); treeseq.tables.num_nodes()],
            left_root: NodeId::NULL,
            above_sample: vec![0; treeseq.tables.num_nodes()],
            left: Position::MIN,
            right: Position::MIN,
            samples: treeseq.samples.as_slice(),
            sample_index_map: vec![NodeId::NULL; treeseq.tables.num_nodes()],
            flags,
            treeseq,
            input_edge_index: 0,
            output_edge_index: 0,
            x: Position(0),
            advanced: false,
        }
    }

    fn init_samples(&mut self) {
        for (i, s) in self.samples.iter().enumerate() {
            if self.sample_index_map[s.0 as usize] != NodeId::NULL {
                panic!("Duplicate samples passed to Tree!");
            }
            self.sample_index_map[s.0 as usize] = NodeId::from(i);
            if let Some(row) = self.topology.get_mut(s.0 as usize) {
                row.left_sample = self.sample_index_map[s.0 as usize];
                row.right_sample = self.sample_index_map[s.0 as usize];
                row.leaf_counts = 1;
                self.above_sample[s.0 as usize] = 1;

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

    fn update_incoming_leaf_count(&mut self, parent: NodeId, child: NodeId) {
        let mut u = parent;
        let lc = self.topology[child.0 as usize].leaf_counts;
        if lc == 0 {
            return;
        }
        while u != NodeId::NULL {
            self.topology[u.0 as usize].leaf_counts += lc;
            u = self.topology[u.0 as usize].parent;
        }
    }

    fn update_incoming_roots(&mut self, parent: NodeId, child: NodeId, lsib: NodeId, rsib: NodeId) {
        if self.above_sample[child.0 as usize] > 0 {
            let mut x = parent;
            let mut root = x;
            let mut above_sample = false;

            while x != NodeId::NULL && !above_sample {
                above_sample = self.above_sample[x.0 as usize] > 0;
                // c is above_sample and p is c's parent.
                // Thus, all parents to p are above_sample, too.
                self.above_sample[x.0 as usize] = 1;
                root = x;
                x = self.topology[x.0 as usize].parent;
            }

            if !above_sample {
                // If we are here, then the above loop terminated
                // by encountering a NULL node, because above_sample[x]
                // must have been 0 for all x. However, because c is
                // above sample, all nodes encountered have been update
                // to be above_sample as well. Thus, the new value of root
                // replaces c in the root list.

                if lsib != NodeId::NULL {
                    self.topology[lsib.0 as usize].right_sib = root;
                }
                if rsib != NodeId::NULL {
                    self.topology[rsib.0 as usize].left_sib = root;
                }
                self.topology[root.0 as usize].left_sib = lsib;
                self.topology[root.0 as usize].right_sib = rsib;
                self.left_root = root;
            } else {
                // If we are here, then we encountered a node
                // ancestral to c where above_sample == 1.
                // Thus, c can no longer be a root.  If the current
                // p is also a c in a later call to this function, then
                // it may also be removed, etc..
                self.left_root = NodeId::NULL;
                if lsib != NodeId::NULL {
                    self.topology[lsib.0 as usize].right_sib = rsib;
                    self.left_root = lsib;
                }
                if rsib != NodeId::NULL {
                    self.topology[rsib.0 as usize].left_sib = lsib;
                    self.left_root = rsib;
                }
            }
        }
    }

    fn update_outgoing_leaf_count(&mut self, parent: NodeId, child: NodeId) {
        let mut u = parent;
        let lc = self.topology[child.0 as usize].leaf_counts;
        if lc == 0 {
            return;
        }
        while u != NodeId::NULL {
            self.topology[u.0 as usize].leaf_counts -= lc;
            u = self.topology[u.0 as usize].parent;
        }
    }

    fn update_outgoing_roots(&mut self, parent: NodeId, child: NodeId) {
        if self.above_sample[child.0 as usize] == 1 {
            let mut x = parent;
            let mut root = x;
            let mut above_sample = false;

            while x != NodeId::NULL && !above_sample {
                above_sample = self.sample_index_map[x.0 as usize] != NodeId::NULL;
                let mut lc = self.topology[x.0 as usize].left_child;
                while lc != NodeId::NULL && !above_sample {
                    above_sample = above_sample || self.above_sample[lc.0 as usize] > 0;
                    lc = self.topology[lc.0 as usize].left_sib;
                }
                if above_sample {
                    self.above_sample[x.0 as usize] = 1;
                }
                root = x;
                x = self.topology[x.0 as usize].parent;
            }

            // Now, root refers to the most ancient
            // ancestor of parent found in the above loop
            if !above_sample {
                // remove root from list of roots
                let lroot = self.topology[root.0 as usize].left_sib;
                let rroot = self.topology[root.0 as usize].right_sib;
                self.left_root = NodeId::NULL;
                if lroot != NodeId::NULL {
                    self.topology[lroot.0 as usize].right_sib = rroot;
                    self.left_root = lroot;
                }
                if rroot != NodeId::NULL {
                    self.topology[rroot.0 as usize].left_sib = lroot;
                    self.left_root = rroot;
                }
                self.topology[root.0 as usize].left_sib = NodeId::NULL;
                self.topology[root.0 as usize].right_sib = NodeId::NULL;
            }
            if self.left_root != NodeId::NULL {
                let lroot = self.topology[self.left_root.0 as usize].left_sib;
                if lroot != NodeId::NULL {
                    self.topology[lroot.0 as usize].right_sib = child;
                }
                self.topology[child.0 as usize].left_sib = lroot;
                self.topology[self.left_root.0 as usize].left_sib = child;
            }
            self.topology[child.0 as usize].right_sib = self.left_root;
            self.left_root = child;
        }
    }

    fn update_samples_list(&mut self, node: NodeId) {
        assert!(self.flags.contains(TreeFlags::TRACK_SAMPLES));

        let sample_map = self.sample_index_map.as_slice();
        let topo = self.topology.as_mut_slice();
        let mut n = node;

        while n != NodeId::NULL {
            let sample_index = sample_map[n.0 as usize];
            if sample_index != NodeId::NULL {
                topo[n.0 as usize].right_sample = topo[n.0 as usize].left_sample;
            } else {
                topo[n.0 as usize].left_sample = NodeId::NULL;
                topo[n.0 as usize].right_sample = NodeId::NULL;
            }

            let mut v = topo[n.0 as usize].left_child;
            while v != NodeId::NULL {
                if topo[v.0 as usize].left_sample != NodeId::NULL {
                    assert!(topo[v.0 as usize].right_sample != NodeId::NULL);
                    if topo[n.0 as usize].left_sample == NodeId::NULL {
                        topo[n.0 as usize].left_sample = topo[v.0 as usize].left_sample;
                    } else {
                        let nright = topo[n.0 as usize].right_sample.0 as usize;
                        let vleft = topo[v.0 as usize].left_sample;
                        topo[nright].next_sample = vleft;
                    }
                    topo[n.0 as usize].right_sample = topo[v.0 as usize].right_sample;
                }
                v = topo[v.0 as usize].right_sib;
            }
            n = topo[n.0 as usize].parent;
        }
    }

    fn id_in_range<N: Into<NodeId>>(&self, u: N) -> TreesResult<()> {
        let n = u.into();
        if n < 0 || (n.0 as usize) >= self.num_nodes() {
            Err(TreesError::NodeIdOutOfRange)
        } else {
            Ok(())
        }
    }

    fn new(treeseq: &'treeseq TreeSequence, flags: TreeFlags) -> Self {
        let mut rv = Self::new_internal(treeseq, flags);
        rv.init_samples();
        rv.left_root = rv.samples[0];
        rv
    }

    /// Return an [`Iterator`] over all nodes in the tree.
    ///
    /// # Parameters
    ///
    /// * `order`: A value from [`NodeTraversalOrder`] specifying the
    ///   iteration order.
    ///
    /// # Errors
    ///
    /// TODO
    // Return value is dyn for later addition of other traversal orders
    // FIXME This should be allowed to error
    pub fn traverse_nodes(
        &self,
        order: NodeTraversalOrder,
    ) -> Box<dyn Iterator<Item = NodeId> + '_> {
        match order {
            NodeTraversalOrder::Preorder => Box::new(PreorderNodeIterator::new(self)),
        }
    }

    /// Return the length of this tree along the genome.
    pub fn span(&self) -> i64 {
        self.right.0 - self.left.0
    }

    /// Return the `[left, right)` [`Position`] for
    /// which this tree is the genealogy.
    pub fn range(&self) -> (Position, Position) {
        (self.left, self.right)
    }

    /// Calculate the total length of the tree via a preorder traversal.
    ///
    /// # Parameters
    ///
    /// * `by_span`: if `true`, multiply the return value by [`Tree::span`].
    pub fn total_branch_length(&self, by_span: bool) -> Result<Time, TreesError> {
        let nt = self.treeseq.tables.nodes_.as_slice();
        let mut b: Time = Time(0.);
        for n in self.traverse_nodes(NodeTraversalOrder::Preorder) {
            let p = self.parent(n)?;
            if p != NodeId::NULL {
                b.0 += nt[n.0 as usize].time.0 - nt[p.0 as usize].time.0;
            }
        }

        match by_span {
            true => Ok(Time(b.0 * (self.span() as f64))),
            false => Ok(b),
        }
    }

    /// Return an [`Iterator`] from the node `u` to the root of the tree,
    /// travering all parent nodes.
    ///
    /// # Errors
    ///
    /// [`TreesError::NodeIdOutOfRange`] if `u` is out of range.
    pub fn parents<N: Into<NodeId> + Copy>(
        &self,
        u: N,
    ) -> Result<impl Iterator<Item = NodeId> + '_, TreesError> {
        self.id_in_range(u)?;
        Ok(ParentsIterator::new(self, u.into()))
    }

    /// Return an [`Iterator`] over the children of node `u`.
    ///
    /// # Errors
    ///
    /// [`TreesError::NodeIdOutOfRange`] if `u` is out of range.
    pub fn children<N: Into<NodeId> + Copy>(
        &self,
        u: N,
    ) -> Result<impl Iterator<Item = NodeId> + '_, TreesError> {
        self.id_in_range(u)?;
        Ok(ChildIterator::new(self, u.into()))
    }

    /// Return an [`Iterator`] over the roots of the tree.
    ///
    /// # Note
    ///
    /// For a tree with multiple roots, the iteration starts
    /// at the left root.
    pub fn roots(&self) -> impl Iterator<Item = NodeId> + '_ {
        RootIterator::new(self)
    }

    /// Return all roots as a vector.
    pub fn roots_to_vec(&self) -> Vec<NodeId> {
        let mut v = vec![];

        for r in self.roots() {
            v.push(r);
        }

        v
    }

    /// Return a slice of the samples in this tree.
    pub fn sample_nodes(&self) -> &[NodeId] {
        self.samples
    }

    /// Return an [`Iterator`] over the sample nodes descending from node `u`.
    ///
    ///
    /// # Note
    ///
    /// If `u` is itself a sample, then it is included in the values returned.
    ///
    /// # Errors
    ///
    /// [`TreesError::NodeIdOutOfRange`] if `u` is out of range.
    ///
    /// [`TreesError::NotTrackingSamples`] if [`TreeFlags::TRACK_SAMPLES`] was not used
    /// to initialize `self`.
    pub fn samples<N: Into<NodeId> + Copy>(
        &self,
        u: N,
    ) -> Result<impl Iterator<Item = NodeId> + '_, TreesError> {
        if !self.flags.contains(TreeFlags::TRACK_SAMPLES) {
            Err(TreesError::NotTrackingSamples)
        } else {
            Ok(SamplesIterator::new(self, u.into()))
        }
    }

    /// The number of nodes in the tree sequence.
    pub fn num_nodes(&self) -> usize {
        assert_eq!(self.topology.len(), self.treeseq.tables.num_nodes());
        self.treeseq.tables.num_nodes()
    }

    /// Return the parent of node `u`.
    pub fn parent<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.parent)
    }

    /// Return the left child of node `u`.
    pub fn left_child<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.left_child)
    }

    /// Return the right child of node `u`.
    pub fn right_child<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.right_child)
    }

    /// Return the left sibling of node `u`.
    pub fn left_sib<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.left_sib)
    }

    /// Return the right sibling of node `u`.
    pub fn right_sib<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.right_sib)
    }

    /// Return the left sample of node `u`.
    pub fn left_sample<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        if !self.flags.contains(TreeFlags::TRACK_SAMPLES) {
            return Err(TreesError::NotTrackingSamples);
        }
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.left_sample)
    }
    //
    /// Return the next sample after node `u`.
    pub fn next_sample<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        if !self.flags.contains(TreeFlags::TRACK_SAMPLES) {
            return Err(TreesError::NotTrackingSamples);
        }
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.next_sample)
    }

    /// Return the right sample of node `u`.
    pub fn right_sample<N: Into<NodeId> + Copy>(&self, u: N) -> TreesResult<NodeId> {
        if !self.flags.contains(TreeFlags::TRACK_SAMPLES) {
            return Err(TreesError::NotTrackingSamples);
        }
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.into().0 as usize) }.right_sample)
    }
}

/// Left-to-right iteration of trees.
impl<'treeseq> streaming_iterator::StreamingIterator for Tree<'treeseq> {
    type Item = Tree<'treeseq>;

    // TODO: if tables are validated when TreeSequence is created,
    // then the accesses below can be unchecked.
    fn advance(&mut self) {
        let tables = &self.treeseq.tables;
        let edge_table = self.treeseq.tables.edges_.as_slice();
        let edge_input_order = tables.edge_input_order.as_slice();
        let edge_output_order = tables.edge_output_order.as_slice();
        if self.input_edge_index < edge_input_order.len() || self.x < tables.genome_length() {
            for edge_index in edge_output_order[self.output_edge_index..].iter() {
                let current_edge = edge_table[*edge_index];
                if current_edge.right != self.x {
                    break;
                }
                let lsib = self.topology[current_edge.child.0 as usize].left_sib;
                let rsib = self.topology[current_edge.child.0 as usize].right_sib;

                if lsib == NodeId::NULL {
                    self.topology[current_edge.parent.0 as usize].left_child = rsib;
                } else {
                    self.topology[lsib.0 as usize].right_sib = rsib;
                }
                if rsib == NodeId::NULL {
                    self.topology[current_edge.parent.0 as usize].right_child = lsib;
                } else {
                    self.topology[rsib.0 as usize].left_sib = lsib;
                }
                let child_topo = &mut self.topology[current_edge.child.0 as usize];
                child_topo.parent = NodeId::NULL;
                child_topo.left_sib = NodeId::NULL;
                child_topo.right_sib = NodeId::NULL;

                self.update_outgoing_leaf_count(current_edge.parent, current_edge.child);

                if self.flags.contains(TreeFlags::TRACK_SAMPLES) {
                    self.update_samples_list(current_edge.parent);
                }
                self.update_outgoing_roots(current_edge.parent, current_edge.child);
                self.output_edge_index += 1;
            }
            for edge_index in edge_input_order[self.input_edge_index..].iter() {
                let current_edge = edge_table[*edge_index];
                if current_edge.left != self.x {
                    break;
                }
                let rchild = self.topology[current_edge.parent.0 as usize].right_child;
                let lsib = self.topology[current_edge.child.0 as usize].left_sib;
                let rsib = self.topology[current_edge.child.0 as usize].right_sib;

                if rchild == NodeId::NULL {
                    self.topology[current_edge.parent.0 as usize].left_child = current_edge.child;
                    self.topology[current_edge.child.0 as usize].left_sib = NodeId::NULL;
                } else {
                    self.topology[rchild.0 as usize].right_sib = current_edge.child;
                    self.topology[current_edge.child.0 as usize].left_sib = rchild;
                }
                self.topology[current_edge.child.0 as usize].right_sib = NodeId::NULL;
                self.topology[current_edge.child.0 as usize].parent = current_edge.parent;
                self.topology[current_edge.parent.0 as usize].right_child = current_edge.child;

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
            if self.left_root != NodeId::NULL {
                while self.topology[self.left_root.0 as usize].left_sib != NodeId::NULL {
                    self.left_root = self.topology[self.left_root.0 as usize].left_sib;
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
            true => Some(self),
            false => None,
        }
    }
}

/// Error type related to [``TreeSequence``] and [``Tree``].
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum TreesError {
    /// Returned by [``TreeSequence::new``].
    #[error("Tables not indexed.")]
    TablesNotIndexed,
    /// Returned when a [`NodeId`] is not
    /// present in a [`Tree`] or [`TreeSequence`].
    #[error("Node ID out of range")]
    NodeIdOutOfRange,
    /// Returned if a tree sequence is
    /// initialized with no samples.
    #[error("No samples found.")]
    NoSamples,
    /// Returned when there are problems with sample lists.
    #[error("Invalid samples.")]
    InvalidSamples,
    /// Returned if sample lists contain duplicate [`NodeId`].
    #[error("Duplicate samples.")]
    DuplicateSamples,
    /// Returned when information about samples describing
    /// from a node is requested, yet
    /// [`TreeFlags::TRACK_SAMPLES`] is not set}.
    #[error("Not tracking samples.")]
    NotTrackingSamples,
}

/// A tree sequence.
pub struct TreeSequence {
    tables: crate::TableCollection,
    samples: Vec<NodeId>,
    num_trees: u32,
}

/// Result type for operations on trees and tree sequences.
pub type TreesResult<T> = Result<T, TreesError>;

bitflags! {
    /// Bit flags modifying the behavior of [`TreeSequence`]
    /// initialization.
    pub struct TreeSequenceFlags: u32 {
        /// Do not validate tables when creating a [`TreeSequence`]
        const NO_TABLE_VALIDATION = 1 << 0;
    }
}

impl TreeSequence {
    fn new_from_tables(tables: crate::TableCollection) -> Result<Self, Box<dyn std::error::Error>> {
        if !tables.is_indexed() {
            return Err(Box::new(crate::TablesError::TablesNotIndexed));
        }
        let mut samples = vec![];
        for (i, n) in tables.nodes_.iter().enumerate() {
            if n.flags & crate::NodeFlags::IS_SAMPLE.bits() > 0 {
                samples.push(NodeId::from(i));
            }
        }
        if samples.is_empty() {
            Err(Box::new(TreesError::NoSamples))
        } else {
            let num_trees = tables.count_trees()?;
            Ok(Self {
                tables,
                samples,
                num_trees,
            })
        }
    }

    /// Create a new tree sequence from a [`TableCollection`](crate::TableCollection).
    ///
    /// The input tables are consumed, owned by the tree sequence.
    ///
    /// By default, the tables will be validated.
    ///
    /// To disable validation, `flags` should contain
    /// [`TreeSequenceFlags::NO_TABLE_VALIDATION`].
    ///
    /// The list of samples will be populated from the [`node flags`](crate::Node::flags).
    /// Any `flag` containing [`IS_SAMPLE`](crate::NodeFlags::IS_SAMPLE) will be
    /// in the list.
    ///
    /// # Errors
    ///
    /// [`TablesNotIndexed`](crate::TablesError::TablesNotIndexed) if
    /// [`build_indexes`](crate::TableCollection::build_indexes) as not been called.
    ///
    /// [`TablesError`](crate::TablesError) if table validation fails.
    pub fn new(
        tables: crate::TableCollection,
        flags: TreeSequenceFlags,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if !tables.is_indexed() {
            return Err(Box::new(crate::TablesError::TablesNotIndexed));
        }
        if !flags.contains(TreeSequenceFlags::NO_TABLE_VALIDATION) {
            tables.validate(crate::TableValidationFlags::empty())?;
        }
        Self::new_from_tables(tables)
    }

    /// Create a new tree sequence from a table collection
    /// and a list of samples.
    ///
    /// Unlike [`TreeSequence::new`], this function ignores node flags and uses the samples
    /// list instead.
    ///
    /// # Error
    ///
    /// [`TreesError`] if the samples list is invalid.
    pub fn new_with_samples(
        tables: crate::TableCollection,
        samples: &[NodeId],
        flags: TreeSequenceFlags,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if !flags.contains(TreeSequenceFlags::NO_TABLE_VALIDATION) {
            tables.validate(crate::TableValidationFlags::empty())?;
        }
        if samples.is_empty() {
            return Err(Box::new(TreesError::NoSamples));
        }
        let mut nodes = vec![0; tables.nodes_.len()];
        for s in samples {
            if *s == NodeId::NULL {
                return Err(Box::new(TreesError::InvalidSamples));
            }
            if nodes[s.0 as usize] != 0 {
                return Err(Box::new(TreesError::DuplicateSamples));
            }
            nodes[s.0 as usize] = 1;
        }
        for n in tables.nodes_.iter() {
            if n.flags | crate::NodeFlags::IS_SAMPLE.bits() > 0 {
                return Err(Box::new(TreesError::InvalidSamples));
            }
        }
        let num_trees = tables.count_trees()?;
        Ok(Self {
            tables,
            samples: samples.to_vec(),
            num_trees,
        })
    }

    /// Move the underlying [`TableCollection`](crate::TableCollection),
    /// consuming `self`.
    pub fn tables(self) -> crate::TableCollection {
        self.tables
    }

    /// Get a clone of the underlying [`TableCollection`](crate::TableCollection).
    pub fn tables_copy(&self) -> crate::TableCollection {
        self.tables.clone()
    }

    /// Return a streaming iterator over all [`Tree`]
    /// objects in the tree sequence.
    pub fn tree_iterator(&self, flags: TreeFlags) -> Tree<'_> {
        Tree::new(self, flags)
    }

    /// The number of sample nodes
    pub fn sample_nodes(&self) -> &[NodeId] {
        &self.samples
    }

    /// The number of trees in the tree sequence
    pub fn num_trees(&self) -> u32 {
        self.num_trees
    }

    /// Simplify the internal [`TableCollection`](crate::TableCollection).
    ///
    /// # Parameters
    ///
    /// * `samples`: An optional slice of [`NodeId`]
    /// * `flags`: flags to modify the simplification behavior
    ///
    /// # Details
    ///
    /// If `samples` is `None`, then all nodes currently marked
    /// as a sample will be used as a samples list.
    /// Calling this functions without a samples list is only
    /// useful if the tree sequence was created with a valid-but-unsimplified
    /// [`TableCollection`](crate::TableCollection).
    ///
    /// # Returns
    ///
    /// A tuple of [`TableCollection`](crate::TableCollection)
    /// and [`SimplificationOutput`](crate::SimplificationOutput).
    ///
    /// # Errors
    ///
    /// [`SimplificationError`](crate::simplification::SimplificationError)
    ///
    pub fn simplify(
        &self,
        samples: Option<&[NodeId]>,
        flags: crate::SimplificationFlags,
    ) -> Result<(crate::TableCollection, crate::SimplificationOutput), crate::SimplificationError>
    {
        let mut tcopy = self.tables.clone();
        let mut si = crate::SamplesInfo::new();
        match samples {
            Some(x) => si.samples = x.to_vec(),
            None => {
                for (i, n) in self.tables.nodes_.iter().enumerate() {
                    if n.flags | crate::NodeFlags::IS_SAMPLE.bits() > 0 {
                        si.samples.push(NodeId::from(i));
                    }
                }
            }
        }
        let mut output = crate::SimplificationOutput::default();
        crate::simplify_tables_without_state(&si, flags, &mut tcopy, &mut output)?;
        Ok((tcopy, output))
    }
}

// FIXME: these tests are weak.
// We'd be better with a light simulator + direct comparison to tskit.
#[cfg(test)]
mod test_trees {
    use super::*;

    #[test]
    fn test_treeseq_creation_and_table_access() {
        let mut tables = crate::TableCollection::new(100).unwrap();
        tables.add_node(0., 0).unwrap();
        tables
            .add_node_with_flags(1., 0, crate::NodeFlags::IS_SAMPLE.bits())
            .unwrap();
        tables.add_edge(0, 1, 0, 1).unwrap();
        tables
            .build_indexes(crate::IndexTablesFlags::empty())
            .unwrap();

        let ts = TreeSequence::new(tables, TreeSequenceFlags::empty()).unwrap();

        let tref = ts.tables();
        assert_eq!(tref.edges().len(), 1);
    }

    #[test]
    fn test_treeseq_creation_and_tree_creation() {
        let mut tables = crate::TableCollection::new(100).unwrap();
        tables.add_edge(0, 1, 0, 1).unwrap();
        tables
            .add_node_with_flags(0., 0, crate::NodeFlags::IS_SAMPLE.bits())
            .unwrap();
        tables.add_node(1., 0).unwrap();
        tables
            .build_indexes(crate::IndexTablesFlags::empty())
            .unwrap();

        let ts = TreeSequence::new(tables, TreeSequenceFlags::empty()).unwrap();
        let _ = ts.tree_iterator(TreeFlags::empty());
    }

    pub fn make_small_table_collection_two_trees() -> crate::TableCollection {
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
        tables.add_node(0., 0).unwrap();
        tables.add_node(1., 0).unwrap();
        tables
            .add_node_with_flags(2., 0, crate::NodeFlags::IS_SAMPLE.bits())
            .unwrap();
        tables
            .add_node_with_flags(2., 0, crate::NodeFlags::IS_SAMPLE.bits())
            .unwrap();
        tables
            .add_node_with_flags(2., 0, crate::NodeFlags::IS_SAMPLE.bits())
            .unwrap();
        tables
            .add_node_with_flags(2., 0, crate::NodeFlags::IS_SAMPLE.bits())
            .unwrap();
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
        assert_eq!(tables.count_trees().unwrap(), 2);
        tables
    }

    // FIXME: this test is UGLY
    // 1. too much repetition.
    // 2. We should be making loads of random toplogies
    //    and comparing to tskit.
    #[test]
    fn test_two_trees() {
        use streaming_iterator::StreamingIterator;
        let tables = make_small_table_collection_two_trees();
        let treeseq = TreeSequence::new(tables, TreeSequenceFlags::empty()).unwrap();
        assert_eq!(treeseq.samples.len(), 4);

        let mut tree_iter = treeseq.tree_iterator(TreeFlags::TRACK_SAMPLES);
        let mut ntrees = 0;
        // NOTE: all of these expected outputs will
        // change if/when we update to a tree structure
        // using virtual roots.
        let expected_number_of_roots = vec![2, 1];
        let mut expected_root_ids = vec![
            vec![NodeId::from(0)],
            vec![NodeId::from(0), NodeId::from(1)],
        ];
        let mut expected_preorders = vec![vec![], vec![]];
        for i in &[0, 2, 3, 1, 4, 5] {
            expected_preorders[1].push(NodeId::from(*i));
        }
        for i in &[0, 3, 1, 4, 5, 2] {
            expected_preorders[0].push(NodeId::from(*i));
        }
        while let Some(tree) = tree_iter.next() {
            println!("{}", tree.left_root.0);
            if ntrees == 0 {
                let mut nodes = vec![0; tree.num_nodes()];
                for c in tree.children(0).unwrap() {
                    nodes[usize::from(c)] = 1;
                }
                assert_eq!(nodes[2], 1);
                assert_eq!(nodes[3], 1);
                for x in &mut nodes {
                    *x = 0;
                }
                for c in tree.children(1).unwrap() {
                    nodes[usize::from(c)] = 1;
                }
                assert_eq!(nodes[4], 1);
                assert_eq!(nodes[5], 1);

                for p in tree.parents(2).unwrap() {
                    nodes[usize::from(p)] = 1;
                }
                assert_eq!(nodes[0], 1);
                for x in &mut nodes {
                    *x = 0;
                }
                for p in tree.parents(5).unwrap() {
                    nodes[p.0 as usize] = 1;
                }
                assert_eq!(nodes[1], 1);
                for x in &mut nodes {
                    *x = 0;
                }
                let roots = tree.roots_to_vec();
                assert_eq!(roots.len(), 2);
                for r in &roots {
                    nodes[usize::from(*r)] = 1;
                }
                for i in &[0, 1] {
                    assert_eq!(nodes[*i as usize], 1);
                }

                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(0).unwrap() {
                    nodes[usize::from(s)] = 1;
                }
                for i in &[2, 3] {
                    assert_eq!(nodes[*i as usize], 1);
                }
                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(1).unwrap() {
                    nodes[usize::from(s)] = 1;
                }
                for i in &[4, 5] {
                    assert_eq!(nodes[*i as usize], 1);
                }
                for x in &mut nodes {
                    *x = 0;
                }
            } else if ntrees == 1 {
                let mut nodes = vec![0; tree.num_nodes()];
                for c in tree.children(0).unwrap() {
                    nodes[usize::from(c)] = 1;
                }
                assert_eq!(nodes[1], 1);
                assert_eq!(nodes[3], 1);
                for x in &mut nodes {
                    *x = 0;
                }
                for c in tree.children(1).unwrap() {
                    nodes[usize::from(c)] = 1;
                }
                assert_eq!(nodes[2], 1);
                assert_eq!(nodes[4], 1);
                assert_eq!(nodes[5], 1);
                for x in &mut nodes {
                    *x = 0;
                }
                let roots = tree.roots_to_vec();
                assert_eq!(roots.len(), 1);
                for r in &roots {
                    nodes[usize::from(*r)] = 1;
                }
                for i in &[0] {
                    assert_eq!(nodes[*i as usize], 1);
                }
                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(0).unwrap() {
                    nodes[usize::from(s)] = 1;
                }
                for s in tree.sample_nodes() {
                    assert_eq!(nodes[usize::from(*s)], 1);
                }
                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(1).unwrap() {
                    nodes[usize::from(s)] = 1;
                }
                for s in &[2, 4, 5] {
                    assert_eq!(nodes[*s as usize], 1);
                }
            }
            let mut num_roots = 0;
            let eroot_ids = expected_root_ids.pop().unwrap();
            for (i, r) in tree.roots().enumerate() {
                num_roots += 1;
                assert_eq!(r, eroot_ids[i]);
            }
            assert_eq!(expected_number_of_roots[ntrees as usize], num_roots);

            let expected_preorder = expected_preorders.pop().unwrap();

            for (i, n) in tree
                .traverse_nodes(NodeTraversalOrder::Preorder)
                .enumerate()
            {
                assert_eq!(expected_preorder[i], n);
            }

            // Check that each sample node contains itself
            // when iterating over samples.
            for s in tree.sample_nodes() {
                for i in tree.samples(*s).unwrap() {
                    assert_eq!(i, *s);
                }
            }
            ntrees += 1;
        }
        assert_eq!(ntrees, 2);
    }
}

#[cfg(test)]
mod test_treeseq_encapsulation {
    use super::*;

    struct MyStruct {
        tables: crate::TableCollection,
    }

    impl MyStruct {
        fn treeseq(&self) -> TreeSequence {
            TreeSequence::new(self.tables.clone(), TreeSequenceFlags::empty()).unwrap()
        }
    }

    #[test]
    fn test_create_treeseq() {
        use streaming_iterator::StreamingIterator;
        let tables = test_trees::make_small_table_collection_two_trees();
        let mystruct = MyStruct { tables };
        let treeseq = mystruct.treeseq();

        let mut tree_iter = treeseq.tree_iterator(TreeFlags::TRACK_SAMPLES);
        let mut ntrees = 0;
        while tree_iter.next().is_some() {
            ntrees += 1;
        }
        assert_eq!(ntrees, 2);
    }
}
