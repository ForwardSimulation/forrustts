use crate::newtypes::{IdType, NodeId, Position, Time};
use bitflags::bitflags;

bitflags! {
    pub struct TreeFlags: u32 {
        const NONE = 0;
        const TRACK_SAMPLES = 1 << 0;
    }
}

#[derive(Copy, Clone)]
pub struct TopologyData {
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
    root_stack: Vec<NodeId>,
    node_stack: Vec<NodeId>,
    tree: &'a Tree<'a>,
    current_node_: Option<NodeId>,
}

impl<'a> PreorderNodeIterator<'a> {
    fn new(tree: &'a Tree) -> Self {
        let mut rv = PreorderNodeIterator {
            root_stack: tree.roots_to_vec(),
            node_stack: vec![],
            tree,
            current_node_: None,
        };
        rv.root_stack.reverse();
        if let Some(root) = rv.root_stack.pop() {
            rv.node_stack.push(root);
        }
        rv
    }
}

impl NodeIterator for PreorderNodeIterator<'_> {
    fn next_node(&mut self) {
        self.current_node_ = self.node_stack.pop();
        match self.current_node_ {
            Some(u) => {
                let mut c = self.tree.left_child(u).unwrap();
                while c != NodeId::NULL {
                    self.node_stack.push(c);
                    c = self.tree.right_sib(c).unwrap();
                }
            }
            None => {
                if let Some(r) = self.root_stack.pop() {
                    self.current_node_ = Some(r);
                }
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

pub struct Tree<'treeseq> {
    topology: Vec<TopologyData>,
    left_root: NodeId,
    above_sample: Vec<i8>,
    left: Position,
    right: Position,
    samples: &'treeseq [NodeId],
    sample_index_map: Vec<NodeId>, // TODO: decide if this is better as usize.
    flags: TreeFlags,
    treeseq: &'treeseq TreeSequence<'treeseq>,
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
            self.sample_index_map[s.0 as usize] = NodeId(i as IdType);
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

    fn id_in_range(&self, u: NodeId) -> TreesResult<()> {
        if u < 0 || (u.0 as usize) >= self.num_nodes() {
            Err(TreesError::NodeIdOutOfRange)
        } else {
            Ok(())
        }
    }

    pub fn new(treeseq: &'treeseq TreeSequence, flags: TreeFlags) -> Self {
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

    pub fn span(&self) -> i64 {
        self.right.0 - self.left.0
    }

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
    pub fn parents(&self, u: NodeId) -> Result<impl Iterator<Item = NodeId> + '_, TreesError> {
        self.id_in_range(u)?;
        Ok(ParentsIterator::new(self, u))
    }

    /// Return an [`Iterator`] over the children of node `u`.
    ///
    /// # Errors
    ///
    /// [`TreesError::NodeIdOutOfRange`] if `u` is out of range.
    pub fn children(&self, u: IdType) -> Result<impl Iterator<Item = IdType> + '_, TreesError> {
        self.id_in_range(u)?;
        Ok(ChildIterator::new(self, u))
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
    pub fn samples(&self, u: NodeId) -> Result<impl Iterator<Item = NodeId> + '_, TreesError> {
        if !self.flags.contains(TreeFlags::TRACK_SAMPLES) {
            Err(TreesError::NotTrackingSamples)
        } else {
            Ok(SamplesIterator::new(self, u))
        }
    }

    pub fn num_nodes(&self) -> usize {
        assert_eq!(self.topology.len(), self.treeseq.tables.num_nodes());
        self.treeseq.tables.num_nodes()
    }

    pub fn parent(&self, u: NodeId) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.0 as usize) }.parent)
    }

    pub fn left_child(&self, u: NodeId) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.0 as usize) }.left_child)
    }

    pub fn right_child(&self, u: NodeId) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.0 as usize) }.right_child)
    }

    pub fn left_sib(&self, u: NodeId) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.0 as usize) }.left_sib)
    }

    pub fn right_sib(&self, u: NodeId) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.0 as usize) }.right_sib)
    }

    pub fn left_sample(&self, u: NodeId) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.0 as usize) }.left_sample)
    }

    pub fn right_sample(&self, u: NodeId) -> TreesResult<NodeId> {
        self.id_in_range(u)?;
        // SAFETY: just checked the range.
        Ok(unsafe { self.topology.get_unchecked(u.0 as usize) }.right_sample)
    }
}

/// Left-to-right iteration of trees.
impl<'treeseq> streaming_iterator::StreamingIterator for Tree<'treeseq> {
    type Item = Tree<'treeseq>;

    // TODO: if tables are validated when TreeSequence is created,
    // then the accesses below can be unchecked.
    fn advance(&mut self) {
        let tables = self.treeseq.tables;
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
    /// Raised by [``TreeSequence::new``].
    #[error("Tables not indexed.")]
    TablesNotIndexed,
    #[error("Node ID out of range")]
    NodeIdOutOfRange,
    #[error("No samples found.")]
    NoSamples,
    #[error("Invalid samples.")]
    InvalidSamples,
    #[error("Duplicate samples.")]
    DuplicateSamples,
    #[error("Not tracking samples.")]
    NotTrackingSamples,
}

pub struct TreeSequence<'a> {
    tables: &'a crate::TableCollection,
    samples: Vec<NodeId>,
    num_trees: u32,
}

/// Result type for operations on trees and tree sequences.
pub type TreesResult<T> = Result<T, TreesError>;

bitflags! {
    pub struct TreeSequenceFlags: u32 {
        /// Do not validate tables when creating a [`TreeSequence`]
        const NO_TABLE_VALIDATION = 1 << 0;
    }
}

impl<'a> TreeSequence<'a> {
    fn new_from_tables(
        tables: &'a crate::TableCollection,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if !tables.is_indexed() {
            return Err(Box::new(crate::TablesError::TablesNotIndexed));
        }
        let mut samples = vec![];
        for (i, n) in tables.nodes_.iter().enumerate() {
            if n.flags & crate::NodeFlags::IS_SAMPLE.bits() > 0 {
                samples.push(NodeId(i as IdType));
            }
        }
        if samples.is_empty() {
            Err(Box::new(TreesError::NoSamples))
        } else {
            Ok(Self {
                tables,
                samples,
                num_trees: tables.count_trees().unwrap(),
            })
        }
    }

    pub fn new(
        tables: &'a crate::TableCollection,
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

    pub fn new_with_samples(
        tables: &'a crate::TableCollection,
        samples: &[NodeId],
        flags: TreeSequenceFlags,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if !flags.contains(TreeSequenceFlags::NO_TABLE_VALIDATION) {
            tables.validate(crate::TableValidationFlags::empty())?;
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
        Ok(Self {
            tables,
            samples: samples.to_vec(),
            num_trees: tables.count_trees().unwrap(),
        })
    }

    pub fn tables(&self) -> &'a crate::TableCollection {
        self.tables
    }

    pub fn tree_iterator(&self, flags: TreeFlags) -> Tree<'_> {
        Tree::new(self, flags)
    }

    pub fn sample_nodes(&self) -> &[NodeId] {
        &self.samples
    }

    pub fn num_trees(&self) -> u32 {
        self.num_trees
    }

    pub fn simplify(
        &self,
        samples: Option<&[NodeId]>,
        flags: crate::SimplificationFlags,
    ) -> Result<(crate::TableCollection, crate::SimplificationOutput), crate::ForrusttsError> {
        let mut tcopy = self.tables.clone();
        let mut si = crate::SamplesInfo::new();
        match samples {
            Some(x) => si.samples = x.to_vec(),
            None => {
                for (i, n) in self.tables.nodes_.iter().enumerate() {
                    if n.flags | crate::NodeFlags::IS_SAMPLE.bits() > 0 {
                        si.samples.push(NodeId(i as IdType));
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

        let ts = TreeSequence::new(&tables, TreeSequenceFlags::empty()).unwrap();

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

        let ts = TreeSequence::new(&tables, TreeSequenceFlags::empty()).unwrap();
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
        let treeseq = TreeSequence::new(&tables, TreeSequenceFlags::empty()).unwrap();
        assert_eq!(treeseq.samples.len(), 4);

        let mut tree_iter = treeseq.tree_iterator(TreeFlags::TRACK_SAMPLES);
        let mut ntrees = 0;
        while let Some(tree) = tree_iter.next() {
            if ntrees == 0 {
                let mut nodes = vec![0; tree.num_nodes()];
                for c in tree.children(0).unwrap() {
                    nodes[c as usize] = 1;
                }
                assert_eq!(nodes[2], 1);
                assert_eq!(nodes[3], 1);
                for x in &mut nodes {
                    *x = 0;
                }
                for c in tree.children(1).unwrap() {
                    nodes[c as usize] = 1;
                }
                assert_eq!(nodes[4], 1);
                assert_eq!(nodes[5], 1);

                for p in tree.parents(2).unwrap() {
                    nodes[p as usize] = 1;
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
                    nodes[r.0 as usize] = 1;
                }
                for i in &[0, 1] {
                    assert_eq!(nodes[*i as usize], 1);
                }

                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(0).unwrap() {
                    nodes[s as usize] = 1;
                }
                for i in &[2, 3] {
                    assert_eq!(nodes[*i as usize], 1);
                }
                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(1).unwrap() {
                    nodes[s as usize] = 1;
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
                    nodes[c as usize] = 1;
                }
                assert_eq!(nodes[1], 1);
                assert_eq!(nodes[3], 1);
                for x in &mut nodes {
                    *x = 0;
                }
                for c in tree.children(1).unwrap() {
                    nodes[c as usize] = 1;
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
                    nodes[r.0 as usize] = 1;
                }
                for i in &[0] {
                    assert_eq!(nodes[*i as usize], 1);
                }
                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(0).unwrap() {
                    nodes[s as usize] = 1;
                }
                for s in tree.sample_nodes() {
                    assert_eq!(nodes[*s as usize], 1);
                }
                for x in &mut nodes {
                    *x = 0;
                }
                for s in tree.samples(1).unwrap() {
                    nodes[s as usize] = 1;
                }
                for s in &[2, 4, 5] {
                    assert_eq!(nodes[*s as usize], 1);
                }
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
        fn treeseq(&self) -> TreeSequence<'_> {
            TreeSequence::new(&self.tables, TreeSequenceFlags::empty()).unwrap()
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
