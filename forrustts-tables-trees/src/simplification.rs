use crate::nested_forward_list::NestedForwardList;
use crate::nested_forward_list::NULL_INDEX;
use crate::newtypes::{NodeId, Position, SiteId, TablesIdInteger, Time};
use crate::tables::*;
use crate::traits::private_traits::TableIdPrivate;
use crate::Segment;
use bitflags::bitflags;
use thiserror::Error;

/// Error type returned by tree sequence
/// simplification
#[derive(Error, Debug, PartialEq)]
pub enum SimplificationError {
    #[error("{0:?}")]
    ErrorMessage(String),
    #[error("{0:?}")]
    ListError(#[from] crate::nested_forward_list::NestedForwardListError),
    #[error("{0:?}")]
    TableValidationError(#[from] crate::TablesError),
}

struct SegmentOverlapper {
    segment_queue: Vec<Segment>,
    overlapping: Vec<Segment>,
    left: Position,
    right: Position,
    qbeg: usize,
    qend: usize,
    num_overlaps: usize,
}

struct MutationNodeMapEntry {
    location: usize,
    output_node: NodeId,
    position: Position,
}

impl SegmentOverlapper {
    fn set_partition(&mut self) -> Position {
        let mut tright = Position::MAX;
        let mut b: usize = 0;

        for i in 0..self.num_overlaps {
            if self.overlapping[i].right > self.left {
                self.overlapping[b] = self.overlapping[i];
                tright = std::cmp::min(tright, self.overlapping[b].right);
                b += 1;
            }
        }

        self.num_overlaps = b;

        tright
    }

    // Public interface below

    const fn new() -> SegmentOverlapper {
        SegmentOverlapper {
            segment_queue: vec![],
            overlapping: vec![],
            left: Position(0),
            right: Position::MAX,
            qbeg: std::usize::MAX,
            qend: std::usize::MAX,
            num_overlaps: std::usize::MAX,
        }
    }

    fn init(&mut self) {
        self.qbeg = 0;
        self.qend = self.segment_queue.len() - 1;
        assert!(self.qend < self.segment_queue.len());
        self.num_overlaps = 0;
        self.overlapping.clear();
    }

    fn enqueue(&mut self, left: Position, right: Position, node: NodeId) {
        self.segment_queue.push(Segment { left, right, node });
    }

    fn finalize_queue(&mut self, maxlen: Position) {
        self.segment_queue.sort_by(|a, b| a.left.cmp(&b.left));
        self.segment_queue.push(Segment {
            left: maxlen,
            right: Position(maxlen.0 + 1),
            node: NodeId::NULL,
        });
    }

    fn advance(&mut self) -> bool {
        let mut rv = false;

        self.left = self.right;
        if self.qbeg < self.qend {
            let mut tright = self.set_partition();
            if self.num_overlaps == 0 {
                self.left = self.segment_queue[self.qbeg].left;
            }
            for seg in self
                .segment_queue
                .iter()
                .skip(self.qbeg)
                .take(self.qend - self.qbeg)
            {
                if seg.left == self.left {
                    tright = std::cmp::min(tright, seg.right);
                    // NOTE: I wonder how efficient this is vs C++?
                    self.overlapping.insert(self.num_overlaps, *seg);
                    self.num_overlaps += 1;
                    self.qbeg += 1;
                } else {
                    break;
                }
            }
            self.right = std::cmp::min(self.segment_queue[self.qbeg].left, tright);
            rv = true;
        } else {
            self.right = Position::MAX;
            let tright = self.set_partition();
            if self.num_overlaps > 0 {
                self.right = tright;
                rv = true
            }
        }

        rv
    }

    fn get_left(&self) -> Position {
        self.left
    }

    fn get_right(&self) -> Position {
        self.right
    }

    fn clear_queue(&mut self) {
        self.segment_queue.clear();
    }
}

type AncestryList = NestedForwardList<Segment>;

// For each input node, we keep a list of locations (usize)
// in the input mutation table and the output id for each node,
// which is initialized to NULL_ID. The position helps us not
// have to possibly refer back to the input site table multiple times.
type MutationNodeMap = NestedForwardList<MutationNodeMapEntry>;

fn prep_mutation_node_map(
    num_nodes: usize,
    mutations: &[MutationRecord],
    sites: &[Site],
    mutation_node_map: &mut MutationNodeMap,
) {
    mutation_node_map.clear();
    mutation_node_map.reset(num_nodes);
    assert_eq!(mutation_node_map.len(), num_nodes);

    for (location, m) in mutations.iter().enumerate() {
        let position = sites[m.site.0 as usize].position;
        mutation_node_map
            .extend(
                m.node.0, // input node id
                MutationNodeMapEntry {
                    location,
                    output_node: NodeId::NULL,
                    position,
                },
            )
            .unwrap();
    }
}

fn record_site(
    position: Position,
    sites: &[Site],
    mutation: &mut MutationRecord,
    new_site_table: &mut SiteTable,
) {
    debug_assert_eq!(position, sites[mutation.site.0 as usize].position);
    if new_site_table.is_empty() || new_site_table[new_site_table.len() - 1].position != position {
        new_site_table.push(sites[mutation.site.0 as usize].clone());
    }

    mutation.site = SiteId((new_site_table.len() - 1) as TablesIdInteger);
}

// This behavior is equivalent to tskit's FILTER_SITES
// option when simplifying.
fn generate_output_site_mutation_tables(
    mutation_node_map: &MutationNodeMap,
    input_sites: &mut SiteTable,
    input_mutations: &mut MutationTable,
    output_sites: &mut SiteTable,
    extinct_mutations: &mut Vec<usize>,
) {
    debug_assert!(output_sites.is_empty());
    debug_assert!(extinct_mutations.is_empty());

    let mut current_input_mutation: usize = 0;

    for site_id in 0..input_sites.len() as TablesIdInteger {
        while current_input_mutation < input_mutations.len() {
            let input_mutation = &input_mutations[current_input_mutation];
            if input_mutation.site != site_id {
                break;
            }

            for val in mutation_node_map.values_iter(input_mutation.node.0) {
                if input_mutations[val.location].site == site_id {
                    if val.output_node != NodeId::NULL {
                        record_site(
                            val.position,
                            input_sites,
                            &mut input_mutations[val.location],
                            output_sites,
                        );
                        input_mutations[val.location].node = val.output_node;
                    } else {
                        input_mutations[val.location].node = NodeId::NULL;
                        extinct_mutations.push(val.location);
                    }
                }
            }

            current_input_mutation += 1;
        }
    }

    std::mem::swap(input_sites, output_sites);
    input_mutations.retain(|m| m.node != NodeId::NULL);
}

fn map_mutation_output_nodes(
    input_id: NodeId,
    output_id: NodeId,
    left: Position,
    right: Position,
    mutation_node_map: &mut MutationNodeMap,
) {
    debug_assert!(
        input_id.0 as usize <= mutation_node_map.len(),
        "{} {}",
        input_id.0,
        mutation_node_map.len()
    );
    let mut list_head = mutation_node_map.head(input_id.0).unwrap();

    while list_head != NULL_INDEX {
        let value = mutation_node_map.fetch_mut(list_head).unwrap();
        if left <= value.position && value.position < right {
            // We cannot assert that .output_node != NULL_ID here:
            // If we are adding input roots in, etc., then we
            // may remap the output ID of an alredy-remapped node.
            value.output_node = output_id;
        }
        list_head = mutation_node_map.next(list_head).unwrap();
    }
}

fn find_parent_child_segment_overlap(
    edges: &[Edge],
    edge_index: usize,
    num_edges: usize,
    maxlen: Position,
    u: NodeId,
    ancestry: &mut AncestryList,
    overlapper: &mut SegmentOverlapper,
) -> Result<usize, SimplificationError> {
    overlapper.clear_queue();

    let mut i = edge_index;

    while i < num_edges && edges[i].parent == u {
        let edge = &edges[i];

        for seg in ancestry.values_iter(edge.child.0) {
            if seg.right > edge.left && edge.right > seg.left {
                overlapper.enqueue(
                    std::cmp::max(seg.left, edge.left),
                    std::cmp::min(seg.right, edge.right),
                    seg.node,
                );
            }
        }
        i += 1;
    }
    overlapper.finalize_queue(maxlen);
    Ok(i)
}

fn add_ancestry(
    input_id: NodeId,
    left: Position,
    right: Position,
    node: NodeId,
    mutation_node_map: &mut MutationNodeMap,
    ancestry: &mut AncestryList,
) -> Result<(), SimplificationError> {
    let head = ancestry.head(input_id.0)?;
    if head == NULL_INDEX {
        let seg = Segment { left, right, node };
        ancestry.extend(input_id.0, seg)?;
    } else {
        let last_idx = ancestry.tail(input_id.0)?;
        if last_idx == NULL_INDEX {
            return Err(SimplificationError::ErrorMessage(
                "last_idx is NULL_ID".to_string(),
            ));
        }
        let last = ancestry.fetch_mut(last_idx)?;
        if last.right == left && last.node == node {
            last.right = right;
        } else {
            let seg = Segment { left, right, node };
            ancestry.extend(input_id.0, seg)?;
        }
    }
    map_mutation_output_nodes(input_id, node, left, right, mutation_node_map);
    Ok(())
}

fn buffer_edge(
    left: Position,
    right: Position,
    parent: NodeId,
    child: NodeId,
    temp_edge_buffer: &mut EdgeTable,
) {
    let i = temp_edge_buffer
        .iter()
        .rposition(|e: &Edge| e.child == child);

    match i {
        None => temp_edge_buffer.push(Edge {
            left,
            right,
            parent,
            child,
        }),
        Some(x) => {
            if temp_edge_buffer[x].right == left {
                temp_edge_buffer[x].right = right;
            } else {
                temp_edge_buffer.push(Edge {
                    left,
                    right,
                    parent,
                    child,
                });
            }
        }
    }
}

fn output_buffered_edges(temp_edge_buffer: &mut EdgeTable, new_edges: &mut EdgeTable) -> usize {
    temp_edge_buffer.sort_by(|a, b| a.child.cmp(&b.child));

    // Need to store size here b/c
    // append drains contents of input!!!
    let rv = temp_edge_buffer.len();
    new_edges.append(temp_edge_buffer);

    rv
}

fn record_node(
    input_nodes: &[Node],
    id: NodeId,
    is_sample: bool,
    output_nodes: &mut NodeTable,
    idmap: &mut [NodeId],
) {
    let mut flags = input_nodes[id.0 as usize].flags;
    flags &= !crate::tables::NodeFlags::IS_SAMPLE.bits();
    if is_sample {
        flags |= crate::tables::NodeFlags::IS_SAMPLE.bits();
    }
    output_nodes.push(Node {
        time: input_nodes[id.0 as usize].time,
        deme: input_nodes[id.0 as usize].deme,
        flags,
    });
    idmap[id.0 as usize] = NodeId((output_nodes.len() - 1) as TablesIdInteger);
}

fn merge_ancestors(
    input_nodes: &[Node],
    maxlen: Position,
    parent_input_id: NodeId,
    state: &mut SimplificationBuffers,
    idmap: &mut [NodeId],
) -> Result<(), SimplificationError> {
    let mut output_id = idmap[parent_input_id.0 as usize];
    let is_sample = output_id != NodeId::NULL;

    if is_sample {
        state.ancestry.nullify_list(parent_input_id.0)?;
    }

    let mut previous_right: Position = Position(0);
    let mut ancestry_node: NodeId;
    state.overlapper.init();
    state.temp_edge_buffer.clear();

    while state.overlapper.advance() {
        if state.overlapper.num_overlaps == 1 {
            ancestry_node = state.overlapper.overlapping[0].node;
            if is_sample {
                buffer_edge(
                    state.overlapper.get_left(),
                    state.overlapper.get_right(),
                    output_id,
                    ancestry_node,
                    &mut state.temp_edge_buffer,
                );
                ancestry_node = output_id;
            }
        } else {
            if output_id == NodeId::NULL {
                record_node(
                    input_nodes,
                    parent_input_id,
                    is_sample,
                    &mut state.new_nodes,
                    idmap,
                );
                output_id = idmap[parent_input_id.0 as usize];
            }
            ancestry_node = output_id;
            for o in state
                .overlapper
                .overlapping
                .iter()
                .take(state.overlapper.num_overlaps)
            {
                buffer_edge(
                    state.overlapper.get_left(),
                    state.overlapper.get_right(),
                    output_id,
                    o.node,
                    &mut state.temp_edge_buffer,
                );
            }
        }
        if is_sample && state.overlapper.get_left() != previous_right {
            add_ancestry(
                parent_input_id,
                previous_right,
                state.overlapper.get_left(),
                output_id,
                &mut state.mutation_node_map,
                &mut state.ancestry,
            )?;
        }
        add_ancestry(
            parent_input_id,
            state.overlapper.get_left(),
            state.overlapper.get_right(),
            ancestry_node,
            &mut state.mutation_node_map,
            &mut state.ancestry,
        )?;
        previous_right = state.overlapper.get_right();
    }
    if is_sample && previous_right != maxlen {
        add_ancestry(
            parent_input_id,
            previous_right,
            maxlen,
            output_id,
            &mut state.mutation_node_map,
            &mut state.ancestry,
        )?;
    }

    if output_id != NodeId::NULL {
        let n = output_buffered_edges(&mut state.temp_edge_buffer, &mut state.new_edges);

        if n == 0 && !is_sample {
            assert!(output_id < state.new_nodes.len() as TablesIdInteger);
            state.new_nodes.truncate(output_id.0 as usize);
            idmap[parent_input_id.0 as usize] = NodeId::NULL;
        }
    }
    Ok(())
}

fn record_sample_nodes(
    samples: &[NodeId],
    tables: &TableCollection,
    new_nodes: &mut NodeTable,
    ancestry: &mut AncestryList,
    mutation_node_map: &mut MutationNodeMap,
    idmap: &mut [NodeId],
) -> Result<(), SimplificationError> {
    for sample in samples.iter() {
        assert!(*sample >= 0);
        // NOTE: the following can be debug_assert?
        if *sample == NodeId::NULL {
            return Err(SimplificationError::ErrorMessage(
                "sample node is NULL_ID".to_string(),
            ));
        }
        if idmap[sample.0 as usize] != NodeId::NULL {
            return Err(SimplificationError::ErrorMessage(
                "invalid sample list!".to_string(),
            ));
        }
        record_node(&tables.nodes_, *sample, true, new_nodes, idmap);

        add_ancestry(
            *sample,
            Position(0),
            tables.genome_length(),
            NodeId((new_nodes.len() - 1) as TablesIdInteger),
            mutation_node_map,
            ancestry,
        )?;
    }
    Ok(())
}

fn validate_tables(
    tables: &TableCollection,
    flags: &SimplificationFlags,
) -> Result<(), SimplificationError> {
    if flags.contains(SimplificationFlags::VALIDATE_EDGES) {
        validate_edge_table(tables.genome_length(), tables.edges(), tables.nodes())?;
    }
    Ok(())
}

fn setup_idmap(nodes: &[Node], idmap: &mut Vec<NodeId>) {
    idmap.resize(nodes.len(), NodeId::NULL);
    idmap.iter_mut().for_each(|x| *x = NodeId::NULL);
}

fn setup_simplification(
    samples: &SamplesInfo,
    tables: &TableCollection,
    flags: SimplificationFlags,
    state: &mut SimplificationBuffers,
    output: &mut SimplificationOutput,
) -> Result<(), SimplificationError> {
    validate_tables(tables, &flags)?;

    // TODO: the SimplificationOutput type could use
    // a "init me" function.
    setup_idmap(&tables.nodes_, &mut output.idmap);
    output.extinct_mutations.clear();

    state.clear();
    state.ancestry.reset(tables.num_nodes());

    prep_mutation_node_map(
        tables.nodes().len(),
        &tables.mutations_,
        &tables.sites_,
        &mut state.mutation_node_map,
    );

    record_sample_nodes(
        &samples.samples,
        tables,
        &mut state.new_nodes,
        &mut state.ancestry,
        &mut state.mutation_node_map,
        &mut output.idmap,
    )?;

    debug_assert_eq!(state.mutation_node_map.len(), tables.nodes().len());

    Ok(())
}

fn process_parent(
    u: NodeId,
    (edge_index, num_edges): (usize, usize),
    tables: &TableCollection,
    state: &mut SimplificationBuffers,
    output: &mut SimplificationOutput,
) -> Result<usize, SimplificationError> {
    let edge_i = find_parent_child_segment_overlap(
        &tables.edges_,
        edge_index,
        num_edges,
        tables.genome_length(),
        u,
        &mut state.ancestry,
        &mut state.overlapper,
    )?;

    merge_ancestors(
        &tables.nodes_,
        tables.genome_length(),
        u,
        state,
        &mut output.idmap,
    )?;
    Ok(edge_i)
}

struct ParentLocation {
    parent: NodeId,
    start: usize,
    stop: usize,
}

// TODO: validate input and return errors.
impl ParentLocation {
    fn new(parent: NodeId, start: usize, stop: usize) -> Self {
        ParentLocation {
            parent,
            start,
            stop,
        }
    }
}

fn find_pre_existing_edges(
    tables: &TableCollection,
    edge_buffer_founder_nodes: &[NodeId],
    edge_buffer: &EdgeBuffer,
) -> Result<Vec<ParentLocation>, SimplificationError> {
    let mut alive_with_new_edges: Vec<i32> = vec![];

    for a in edge_buffer_founder_nodes {
        if edge_buffer.0.head(a.0)? != NULL_INDEX {
            alive_with_new_edges.push(a.0);
        }
    }
    if alive_with_new_edges.is_empty() {
        return Ok(vec![]);
    }

    let mut starts = vec![usize::MAX; tables.num_nodes()];
    let mut stops = vec![usize::MAX; tables.num_nodes()];

    for (i, e) in tables.enumerate_edges() {
        if starts[e.parent.0 as usize] == usize::MAX {
            starts[e.parent.0 as usize] = i;
        }
        stops[e.parent.0 as usize] = i + 1;
    }

    let mut rv = vec![];
    for a in alive_with_new_edges {
        rv.push(ParentLocation::new(
            NodeId::new(a),
            starts[a as usize],
            stops[a as usize],
        ));
    }

    rv.sort_by(|a, b| {
        let ta = tables.nodes_[a.parent.0 as usize].time;
        let tb = tables.nodes_[b.parent.0 as usize].time;
        match ta.partial_cmp(&tb) {
            Some(std::cmp::Ordering::Equal) => {
                if a.start == b.start {
                    a.parent.cmp(&b.parent)
                } else {
                    a.start.cmp(&b.start)
                }
            }
            Some(x) => x.reverse(),
            None => panic!("invalid node times"),
        }
    });

    // TODO: this could eventually be called in a debug_assert
    if !rv.is_empty() {
        for (i, _) in rv.iter().enumerate().skip(1) {
            let t0 = tables.nodes_[rv[i - 1].parent.0 as usize].time;
            let t1 = tables.nodes_[rv[i].parent.0 as usize].time;
            if t0 < t1 {
                return Err(SimplificationError::ErrorMessage(
                    "existing edges not properly sorted by time".to_string(),
                ));
            }
        }
    }
    Ok(rv)
}

fn queue_children(
    child: NodeId,
    left: Position,
    right: Position,
    ancestry: &mut AncestryList,
    overlapper: &mut SegmentOverlapper,
) -> Result<(), SimplificationError> {
    for seg in ancestry.values_iter(child.0) {
        if seg.right > left && right > seg.left {
            overlapper.enqueue(
                std::cmp::max(seg.left, left),
                std::cmp::min(seg.right, right),
                seg.node,
            );
        }
    }
    Ok(())
}

fn process_births_from_buffer(
    head: NodeId,
    edge_buffer: &EdgeBuffer,
    state: &mut SimplificationBuffers,
) -> Result<(), SimplificationError> {
    // Have to take references here to
    // make the borrow checker happy.
    let a = &mut state.ancestry;
    let o = &mut state.overlapper;
    for seg in edge_buffer.0.values_iter(head.0) {
        queue_children(seg.node, seg.left, seg.right, a, o).unwrap();
    }
    Ok(())
}

bitflags! {
    /// Boolean flags affecting simplification
    /// behavior.
    ///
    /// # Example
    ///
    /// ```
    /// let e = forrustts_tables_trees::SimplificationFlags::empty();
    /// assert_eq!(e.bits(), 0);
    /// ```
    #[derive(Default)]
    pub struct SimplificationFlags: u32 {
        /// Validate that input edges are sorted
        const VALIDATE_EDGES = 1 << 0;
        /// Validate that input mutations are sorted
        const VALIDATE_MUTATIONS = 1 << 1;
        /// Validate all tables.
        const VALIDATE_ALL = Self::VALIDATE_EDGES.bits | Self::VALIDATE_MUTATIONS.bits;
    }
}

/// Information about samples used for
/// table simpilfication.
#[derive(Default)]
pub struct SamplesInfo {
    /// A list of sample IDs.
    /// Can include both "alive" and
    /// "ancient/remembered/preserved" sample
    /// nodes.
    pub samples: Vec<NodeId>,
    /// When using [``EdgeBuffer``](type.EdgeBuffer.html)
    /// to record transmission
    /// events, this list must contain a list of all node IDs
    /// alive the last time simplification happened. Here,
    /// "alive" means "could leave more descendants".
    /// At the *start* of a simulation, this  should be filled
    /// with a list of "founder" node IDs.
    pub edge_buffer_founder_nodes: Vec<NodeId>,
}

impl SamplesInfo {
    /// Generate a new instance.
    pub fn new() -> Self {
        SamplesInfo {
            samples: vec![],
            edge_buffer_founder_nodes: vec![],
        }
    }
}

/// Useful information output by table
/// simplification.
pub struct SimplificationOutput {
    /// Maps input node ID to output ID.
    /// Values are set to [``NodeId::NULL``](crate::NodeId::NULL)
    /// for input nodes that "simplify out".
    pub idmap: Vec<NodeId>,
    /// Indexes of mutations that simplified out
    pub extinct_mutations: Vec<usize>,
}

impl SimplificationOutput {
    /// Create a new instance.
    pub fn new() -> Self {
        SimplificationOutput {
            idmap: vec![],
            extinct_mutations: vec![],
        }
    }
}

impl Default for SimplificationOutput {
    fn default() -> Self {
        SimplificationOutput::new()
    }
}

/// Holds internal memory used by
/// simplification machinery.
///
/// During simplification, several large
/// memory blocks are required. This type
/// allows those allocations to be re-used
/// in subsequent calls to [``simplify_tables``].
/// Doing so typically improves run times at
/// the cost of higher peak memory consumption.
pub struct SimplificationBuffers {
    new_edges: EdgeTable,
    temp_edge_buffer: EdgeTable,
    new_nodes: NodeTable,
    overlapper: SegmentOverlapper,
    ancestry: AncestryList,
    mutation_map: Vec<MutationNodeMapEntry>,
    mutation_node_map: MutationNodeMap,
    new_site_table: Vec<Site>,
}

impl SimplificationBuffers {
    /// Create a new instance.
    pub const fn new() -> SimplificationBuffers {
        SimplificationBuffers {
            new_edges: EdgeTable::new(),
            temp_edge_buffer: EdgeTable::new(),
            new_nodes: NodeTable::new(),
            overlapper: SegmentOverlapper::new(),
            ancestry: AncestryList::new(),
            mutation_map: vec![],
            mutation_node_map: MutationNodeMap::new(),
            new_site_table: vec![],
        }
    }

    // NOTE: should this be fully pub?
    fn clear(&mut self) {
        self.new_edges.clear();
        self.temp_edge_buffer.clear();
        self.new_nodes.clear();
        self.mutation_map.clear();
        self.mutation_node_map.clear();
        self.new_site_table.clear();
    }
}

impl Default for SimplificationBuffers {
    fn default() -> Self {
        Self::new()
    }
}

/// Simplify a [``TableCollection``].
///
/// # Parameters
///
/// * `samples`:
/// * `flags`: modify the behavior of the simplification algorithm.
/// * `tables`: a [``TableCollection``] to simplify.
/// * `output`: Where simplification output gets written.
///             See [``SimplificationOutput``].
///
/// # Notes
///
/// The input tables must be sorted.
/// See [``TableCollection::sort_tables_for_simplification``].
///
/// It is common to simplify many times during a simulation.
/// To avoid making big allocations each time, see
/// [``simplify_tables``] to keep memory allocations
/// persistent between simplifications.
pub fn simplify_tables_without_state(
    samples: &SamplesInfo,
    flags: SimplificationFlags,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) -> Result<(), SimplificationError> {
    let mut state = SimplificationBuffers::new();
    simplify_tables(samples, flags, &mut state, tables, output)
}

/// Simplify a [``TableCollection``].
///
/// This differs from [``simplify_tables_without_state``] in that the big memory
/// allocations made during simplification are preserved in
/// an instance of [``SimplificationBuffers``].
///
/// # Parameters
///
/// * `samples`:
/// * `flags`: modify the behavior of the simplification algorithm.
/// * `state`: These are the internal data structures used
///            by the simpilfication algorithm.
/// * `tables`: a [``TableCollection``] to simplify.
/// * `output`: Where simplification output gets written.
///             See [``SimplificationOutput``].
///
/// # Notes
///
/// The input tables must be sorted.
/// See [``TableCollection::sort_tables_for_simplification``].
pub fn simplify_tables(
    samples: &SamplesInfo,
    flags: SimplificationFlags,
    state: &mut SimplificationBuffers,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) -> Result<(), SimplificationError> {
    setup_simplification(samples, tables, flags, state, output)?;

    let mut edge_i = 0;
    let num_edges = tables.num_edges();
    let mut new_edges_inserted: usize = 0;
    while edge_i < num_edges {
        edge_i = process_parent(
            tables.edges_[edge_i].parent,
            (edge_i, num_edges),
            tables,
            state,
            output,
        )?;

        if state.new_edges.len() >= 1024 && new_edges_inserted + state.new_edges.len() < edge_i {
            for i in state.new_edges.drain(..) {
                tables.edges_[new_edges_inserted] = i;
                new_edges_inserted += 1;
            }
            assert_eq!(state.new_edges.len(), 0);
        }
    }

    tables.edges_.truncate(new_edges_inserted);
    tables.edges_.append(&mut state.new_edges);
    std::mem::swap(&mut tables.nodes_, &mut state.new_nodes);

    generate_output_site_mutation_tables(
        &state.mutation_node_map,
        &mut tables.sites_,
        &mut tables.mutations_,
        &mut state.new_site_table,
        &mut output.extinct_mutations,
    );

    Ok(())
}

/// Data type used for edge buffering.
/// Simplification of simulated data happens
/// via [``crate::simplify_from_edge_buffer()``].
///
/// # Overview
///
/// The typical tree sequence recording workflow
/// goes like this.  When a new node is "born",
/// we:
///
/// 1. Add a new node to the node table.  This
///    new node is a "child".
/// 2. Add edges to the edge table representing
///    the genomic intervals passed on from
///    various parents to this child node.
///
/// We repeat `1` and `2` for a while, then
/// we [``sort the tables``](crate::TableCollection::sort_tables_for_simplification).
/// After sorting, we [``simplify``](crate::simplify_tables()) the tables.
///
/// We can avoid the sorting step using this type.
/// To start, we record the list of currently-alive nodes
/// [``here``](crate::SamplesInfo::edge_buffer_founder_nodes).
///
/// In a simulation, offspring are generated by birth order.
/// Further, a well-behaved forward simulation is capable of calculating
/// edges from left-to-right along a genome.  Thus, we can
/// [``record``](EdgeBuffer::record_edge) the data representing
/// transmission events.
///
/// After recording for a while, we call
/// [``simplify_from_edge_buffer``](crate::simplify_from_edge_buffer()) to simplify
/// the tables.  After simplification, the client code must re-populate
/// [``the list``](crate::SamplesInfo::edge_buffer_founder_nodes) of alive nodes.
/// Once that is done, we can keep recording, etc..
///
/// # Example
///
/// For a full example of use in simulation,
/// see the source code for `examples/forward_simulation.rs`
///
#[repr(transparent)]
pub struct EdgeBuffer(pub(crate) NestedForwardList<Segment>);

impl EdgeBuffer {
    /// Create a new EdgeBuffer
    pub fn new() -> Self {
        Self(NestedForwardList::new())
    }

    /// Add an edge to the buffer
    pub fn record_edge<
        P: Into<crate::newtypes::NodeId>,
        C: Into<crate::newtypes::NodeId>,
        L: Into<crate::newtypes::Position>,
        R: Into<crate::newtypes::Position>,
    >(
        &mut self,
        parent: P,
        child: C,
        left: L,
        right: R,
    ) -> crate::nested_forward_list::Result<()> {
        self.0.extend(
            parent.into().0,
            crate::segment::Segment::new(left.into(), right.into(), child.into()),
        )?;
        Ok(())
    }
}

impl Default for EdgeBuffer {
    fn default() -> Self {
        Self::new()
    }
}

/// Simplify a [``TableCollection``] from an [``EdgeBuffer``].
///
/// See [``EdgeBuffer``] for discussion.
///
/// # Parameters
///
/// * `samples`: Instance of [``SamplesInfo``]. The field
///              [``SamplesInfo::edge_buffer_founder_nodes``]
///              must be populated. See [``EdgeBuffer``] for details.
/// * `flags`: modify the behavior of the simplification algorithm.
/// * `state`: These are the internal data structures used
///            by the simpilfication algorithm.
/// * `edge_buffer`: An [``EdgeBuffer``] recording births since the last
///                  simplification.
/// * `tables`: a [``TableCollection``] to simplify.
/// * `output`: Where simplification output gets written.
///             See [``SimplificationOutput``].
///
/// # Notes
///
/// The input tables must be sorted.
/// See [``TableCollection::sort_tables_for_simplification``].
///
/// # Limitations
///
/// The simplification code does not currently validate
/// that "buffered" edges do indeed represent a valid sort order.
pub fn simplify_from_edge_buffer(
    samples: &SamplesInfo,
    flags: SimplificationFlags,
    state: &mut SimplificationBuffers,
    edge_buffer: &mut EdgeBuffer,
    tables: &mut TableCollection,
    output: &mut SimplificationOutput,
) -> Result<(), SimplificationError> {
    setup_simplification(samples, tables, flags, state, output)?;

    // Process all edges since the last simplification.
    let mut max_time = Time::MIN;
    for n in samples.edge_buffer_founder_nodes.iter() {
        let nt = tables.node(*n).time;
        max_time = match max_time.partial_cmp(&nt) {
            Some(std::cmp::Ordering::Less) => nt,
            Some(_) => max_time,
            None => panic!("invalid time comparsion"),
        };
    }

    for head in edge_buffer.0.index_rev() {
        let ptime = tables.node(NodeId(head)).time;
        if ptime > max_time
        // Then this is a parent who is:
        // 1. Born since the last simplification.
        // 2. Left offspring
        {
            state.overlapper.clear_queue();
            process_births_from_buffer(NodeId(head), edge_buffer, state)?;
            state.overlapper.finalize_queue(tables.genome_length());
            merge_ancestors(
                &tables.nodes_,
                tables.genome_length(),
                NodeId(head),
                state,
                &mut output.idmap,
            )?;
        } else if ptime <= max_time {
            break;
        }
    }

    let existing_edges =
        find_pre_existing_edges(tables, &samples.edge_buffer_founder_nodes, edge_buffer)?;

    let mut edge_i = 0;
    let num_edges = tables.num_edges();

    for ex in existing_edges {
        while edge_i < num_edges
            && tables.nodes_[tables.edges_[edge_i].parent.0 as usize].time
                > tables.nodes_[ex.parent.0 as usize].time
        {
            edge_i = process_parent(
                tables.edges_[edge_i].parent,
                (edge_i, num_edges),
                tables,
                state,
                output,
            )?;
        }
        if ex.start != usize::MAX {
            while (edge_i as usize) < ex.start
                && tables.nodes_[tables.edges_[edge_i].parent.0 as usize].time
                    >= tables.nodes_[ex.parent.0 as usize].time
            {
                edge_i = process_parent(
                    tables.edges_[edge_i].parent,
                    (edge_i, num_edges),
                    tables,
                    state,
                    output,
                )?;
            }
        }
        // now, handle ex.parent
        state.overlapper.clear_queue();
        if ex.start != usize::MAX {
            while edge_i < ex.stop {
                // TODO: a debug assert or regular assert?
                if tables.edges_[edge_i].parent != ex.parent {
                    return Err(SimplificationError::ErrorMessage(
                        "Unexpected parent node".to_string(),
                    ));
                }
                let a = &mut state.ancestry;
                let o = &mut state.overlapper;
                queue_children(
                    tables.edges_[edge_i].child,
                    tables.edges_[edge_i].left,
                    tables.edges_[edge_i].right,
                    a,
                    o,
                )?;
                edge_i += 1;
            }
            if edge_i < num_edges && tables.edges_[edge_i].parent == ex.parent {
                return Err(SimplificationError::ErrorMessage(
                    "error traversing pre-existing edges for parent".to_string(),
                ));
            }
        }
        process_births_from_buffer(ex.parent, edge_buffer, state)?;
        state.overlapper.finalize_queue(tables.genome_length());
        merge_ancestors(
            &tables.nodes_,
            tables.genome_length(),
            ex.parent,
            state,
            &mut output.idmap,
        )?;
    }

    // Handle remaining edges.
    while edge_i < num_edges {
        edge_i = process_parent(
            tables.edges_[edge_i].parent,
            (edge_i, num_edges),
            tables,
            state,
            output,
        )?;
    }

    std::mem::swap(&mut tables.edges_, &mut state.new_edges);
    std::mem::swap(&mut tables.nodes_, &mut state.new_nodes);
    edge_buffer.0.reset(tables.num_nodes());

    generate_output_site_mutation_tables(
        &state.mutation_node_map,
        &mut tables.sites_,
        &mut tables.mutations_,
        &mut state.new_site_table,
        &mut output.extinct_mutations,
    );

    Ok(())
}

#[cfg(test)]
mod test_samples_info {
    use super::SamplesInfo;

    #[test]
    fn test_default() {
        let s: SamplesInfo = Default::default();
        assert!(s.samples.is_empty());
        assert!(s.edge_buffer_founder_nodes.is_empty());
    }
}

#[cfg(test)]
mod test_simplification_output {
    use super::SimplificationOutput;

    #[test]
    fn test_defaul() {
        let x: SimplificationOutput = Default::default();
        assert!(x.idmap.is_empty());
    }
}

#[cfg(test)]
mod test_simplification_flags {
    use super::SimplificationFlags;

    #[test]
    fn test_empty_simplification_flags() {
        let e = SimplificationFlags::empty();
        assert_eq!(e.bits(), 0);
    }
}

#[cfg(test)]
mod test_simpify_tables {
    use super::*;

    // TODO: we need lots more tests of these validations!

    #[test]
    fn test_simplify_tables_unsorted_edges() {
        let mut tables = TableCollection::new(1000).unwrap();

        tables.add_node(0., 0).unwrap(); // parent
        tables.add_node(1., 0).unwrap(); // child
        tables.add_edge(100, tables.genome_length(), 0, 1).unwrap();
        tables.add_edge(0, 100, 0, 1).unwrap();

        let mut output = SimplificationOutput::new();

        let mut samples = SamplesInfo::new();
        samples.samples.push(1.into());

        let _ = simplify_tables_without_state(
            &samples,
            SimplificationFlags::VALIDATE_ALL,
            &mut tables,
            &mut output,
        )
        .map_or_else(
            |x: SimplificationError| {
                assert_eq!(
                    x,
                    SimplificationError::TableValidationError(TablesError::EdgesNotSortedByLeft),
                )
            },
            |_| panic!(),
        );
    }
}