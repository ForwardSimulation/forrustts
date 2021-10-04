use bitflags::bitflags;
use forrustts::*;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Exp, Uniform};
use std::sync::{Arc, Mutex};
use std::thread;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SimulationError {
    #[error("{0:?}")]
    ErrorMessage(String),
}

// Some of the material below seems like a candidate for a public API,
// but we need to decide here if this package should provide that.
// If so, then many of these types should not be here, as they have nothing
// to do with Wright-Fisher itself, and are instead more general.

#[repr(transparent)]
pub(crate) struct Rng(pub(crate) rgsl::Rng);

impl Rng {
    /// Create a new [`Rng`] with a seed.
    pub fn new(seed: usize) -> Self {
        let mut rng = rgsl::rng::Rng::new(rgsl::rng::algorithms::mt19937()).unwrap();
        rng.set(seed);

        Self(rng)
    }
}

#[derive(Debug, Eq, PartialEq)]
struct Segment {
    left: Position,
    right: Position,
}

impl Segment {
    fn new(left: Position, right: Position) -> Self {
        Self { left, right }
    }
}

impl std::fmt::Display for Segment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Segment({},{})",
            PositionLLType::from(self.left),
            PositionLLType::from(self.right)
        )
    }
}

struct BreakpointIterator<'bp> {
    breakpoints: &'bp [Position],
    left: Option<Position>,
    sequence_length: Position,
    index: usize,
}

enum NextBreakpointIndex {
    Swap,
    Index(usize),
    Done,
}

impl<'bp> BreakpointIterator<'bp> {
    fn new(breakpoints: &'bp [Position], sequence_length: Position) -> Self {
        Self {
            breakpoints,
            left: None,
            sequence_length,
            index: 0,
        }
    }

    fn next_breakpoint(&mut self) -> NextBreakpointIndex {
        let mut t = self.index;

        loop {
            if t < self.breakpoints.len() {
                let candidate = self.breakpoints[t];
                match self.breakpoints[t..].iter().position(|y| *y != candidate) {
                    Some(distance) => {
                        if distance % 2 != 0 {
                            self.index = t + distance;
                            if candidate == 0 {
                                return NextBreakpointIndex::Swap;
                            }
                            return NextBreakpointIndex::Index(t);
                        }
                        t += distance;
                        self.index = t;
                        // Crappy corner case...
                        if candidate == 0 {
                            self.left = Some(candidate);
                            return self.next_breakpoint();
                        }
                    }
                    None => {
                        if self.left.is_some() {
                            self.index = self.breakpoints.len();
                            return NextBreakpointIndex::Done;
                        } else {
                            assert_eq!(self.index, self.breakpoints.len() - 1);
                            self.index = self.breakpoints.len();
                            return NextBreakpointIndex::Index(self.index - 1);
                        }
                    }
                };
            } else {
                return NextBreakpointIndex::Done;
            }
        }
    }
}

// At the start of processing:
// If breakpoints[0] == 0, and it
// occurs an odd number of times,
// Swap.
// If it occurs an even number of times,
// we process from 0, first odd.
// Same if it does not occur.
// We may need to generate this above,
// in next_index_for_breakpoint
#[derive(Debug, Eq, PartialEq)]
enum SegmentAction {
    Swap,
    Process(Segment),
}

impl std::fmt::Display for SegmentAction {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Swap => write!(f, "SegmentAction::Swap"),
            Self::Process(segment) => write!(f, "SegmentAction::Process({})", segment),
        }
    }
}

impl<'bp> Iterator for BreakpointIterator<'bp> {
    type Item = SegmentAction;

    fn next(&mut self) -> Option<SegmentAction> {
        match self.next_breakpoint() {
            NextBreakpointIndex::Swap => Some(SegmentAction::Swap),
            NextBreakpointIndex::Index(index) => match self.left {
                Some(value) => {
                    self.left = Some(self.breakpoints[index]);
                    Some(SegmentAction::Process(Segment::new(
                        value,
                        self.breakpoints[index],
                    )))
                }
                None => {
                    // TODO: is this the best logic here?
                    if self.index < self.breakpoints.len() {
                        self.left = Some(self.breakpoints[index]);
                    }
                    Some(SegmentAction::Process(Segment::new(
                        0.into(),
                        self.breakpoints[index],
                    )))
                }
            },
            NextBreakpointIndex::Done => match self.left {
                Some(left) => {
                    self.left = None;
                    Some(SegmentAction::Process(Segment::new(
                        left,
                        self.sequence_length,
                    )))
                }
                None => None,
            },
        }
    }
}

#[derive(Copy, Clone)]
struct Parent {
    index: usize,
    node0: NodeId,
    node1: NodeId,
}

struct Birth {
    index: usize,
    parent0: Parent,
    parent1: Parent,
}

type VecParent = Vec<Parent>;
type VecBirth = Vec<Birth>;

struct PopulationState {
    pub parents: VecParent,
    pub births: VecBirth,
    pub edge_buffer: EdgeBuffer,
    pub tables: TableCollection,
}

impl PopulationState {
    pub fn new(genome_length: Position) -> Self {
        PopulationState {
            parents: vec![],
            births: vec![],
            edge_buffer: EdgeBuffer::new(),
            tables: TableCollection::new(genome_length).unwrap(),
        }
    }
}

fn deaths_and_parents(psurvival: f64, rng: &mut Rng, pop: &mut PopulationState) {
    pop.births.clear();
    for i in 0..pop.parents.len() {
        let x = rng.0.uniform();
        match x.partial_cmp(&psurvival) {
            Some(std::cmp::Ordering::Greater) => {
                let random_index = rng.0.flat(0.0, pop.parents.len() as f64) as usize;
                let parent0 = pop.parents[random_index];
                let random_index = rng.0.flat(0.0, pop.parents.len() as f64) as usize;
                let parent1 = pop.parents[random_index];
                pop.births.push(Birth {
                    index: i,
                    parent0,
                    parent1,
                });
            }
            Some(_) => (),
            None => (),
        }
    }
}

fn mendel(pnodes: &mut (NodeId, NodeId), rng: &mut Rng) {
    let x: f64 = rng.0.uniform();
    match x.partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Less) => {
            std::mem::swap(&mut pnodes.0, &mut pnodes.1);
        }
        Some(_) => (),
        None => panic!("Unexpected None"),
    }
}

fn crossover_breakpoints(
    recrate: Option<f64>,
    genome_length: Position,
    rng: &mut Rng,
) -> Option<Vec<Position>> {
    match recrate {
        Some(x) => {
            let n = rng.0.poisson(x);
            match n > 0 {
                true => {
                    let mut rv = vec![];
                    for _ in 0..n {
                        rv.push(Position::from(
                            rng.0.flat(0.0, i64::from(genome_length) as f64) as PositionLLType,
                        ));
                    }
                    rv.sort_unstable();
                    rv.push(genome_length); // Sentinel value
                    Some(rv)
                }
                false => None,
            }
        }
        None => None,
    }
}

fn crossover_and_record_edges(
    parent: Parent,
    child: NodeId,
    recrate: Option<f64>,
    recorder: &impl Fn(NodeId, NodeId, (Position, Position), &mut TableCollection, &mut EdgeBuffer),
    rng: &mut Rng,
    tables: &mut TableCollection,
    edge_buffer: &mut EdgeBuffer,
) {
    let mut pnodes = (parent.node0, parent.node1);
    mendel(&mut pnodes, rng);

    match crossover_breakpoints(recrate, tables.genome_length(), rng) {
        Some(x) => {
            let bpiter = BreakpointIterator::new(&x, tables.genome_length());
            for action in bpiter {
                match action {
                    SegmentAction::Swap => (),
                    SegmentAction::Process(segment) => {
                        recorder(
                            pnodes.0,
                            child,
                            (segment.left, segment.right),
                            tables,
                            edge_buffer,
                        );
                    }
                }
                std::mem::swap(&mut pnodes.0, &mut pnodes.1);
            }
        }
        None => {
            recorder(
                pnodes.0,
                child,
                (0.into(), tables.genome_length()),
                tables,
                edge_buffer,
            );
        }
    }
}

fn generate_births(
    recrate: Option<f64>,
    birth_time: Time,
    rng: &mut Rng,
    pop: &mut PopulationState,
    recorder: &impl Fn(NodeId, NodeId, (Position, Position), &mut TableCollection, &mut EdgeBuffer),
) {
    for b in &pop.births {
        // Record 2 new nodes
        let new_node_0: NodeId = pop.tables.add_node(birth_time, 0).unwrap();
        let new_node_1: NodeId = pop.tables.add_node(birth_time, 0).unwrap();

        crossover_and_record_edges(
            b.parent0,
            new_node_0,
            recrate,
            recorder,
            rng,
            &mut pop.tables,
            &mut pop.edge_buffer,
        );
        crossover_and_record_edges(
            b.parent1,
            new_node_1,
            recrate,
            recorder,
            rng,
            &mut pop.tables,
            &mut pop.edge_buffer,
        );

        pop.parents[b.index].index = b.index;
        pop.parents[b.index].node0 = new_node_0;
        pop.parents[b.index].node1 = new_node_1;
    }
}

fn buffer_edges(
    parent: NodeId,
    child: NodeId,
    span: (Position, Position),
    _: &mut TableCollection,
    buffer: &mut EdgeBuffer,
) {
    buffer.record_edge(parent, child, span.0, span.1).unwrap();
}

fn record_edges(
    parent: NodeId,
    child: NodeId,
    span: (Position, Position),
    tables: &mut TableCollection,
    _: &mut EdgeBuffer,
) {
    tables.add_edge(span.0, span.1, parent, child).unwrap();
}

fn fill_samples(parents: &[Parent], samples: &mut SamplesInfo) {
    samples.samples.clear();
    for p in parents {
        samples.samples.push(p.node0);
        samples.samples.push(p.node1);
    }
}

fn sort_and_simplify(
    flags: SimulationFlags,
    simplification_flags: SimplificationFlags,
    samples: &SamplesInfo,
    state: &mut SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut SimplificationOutput,
) {
    if !flags.contains(SimulationFlags::BUFFER_EDGES) {
        pop.tables.sort_tables(TableSortingFlags::empty());
        if flags.contains(SimulationFlags::USE_STATE) {
            simplify_tables(
                samples,
                simplification_flags,
                state,
                &mut pop.tables,
                output,
            )
            .unwrap();
        } else {
            simplify_tables_without_state(samples, simplification_flags, &mut pop.tables, output)
                .unwrap();
        }
    } else {
        simplify_from_edge_buffer(
            samples,
            simplification_flags,
            state,
            &mut pop.edge_buffer,
            &mut pop.tables,
            output,
        )
        .unwrap();
    }
}

fn simplify_and_remap_nodes(
    flags: SimulationFlags,
    simplification_flags: SimplificationFlags,
    samples: &mut SamplesInfo,
    state: &mut SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut SimplificationOutput,
) {
    fill_samples(&pop.parents, samples);
    sort_and_simplify(flags, simplification_flags, samples, state, pop, output);

    for p in &mut pop.parents {
        p.node0 = output.idmap[usize::from(p.node0)];
        p.node1 = output.idmap[usize::from(p.node1)];
        assert!(pop.tables.node(p.node0).flags & NodeFlags::IS_SAMPLE.bits() > 0);
    }

    if flags.contains(SimulationFlags::BUFFER_EDGES) {
        samples.edge_buffer_founder_nodes.clear();
        for p in &pop.parents {
            samples.edge_buffer_founder_nodes.push(p.node0);
            samples.edge_buffer_founder_nodes.push(p.node1);
        }
    }
}

fn validate_simplification_interval(x: i64) -> i64 {
    if x < 1 {
        panic!("simplification_interval must be None or >= 1");
    }
    x
}

bitflags! {
    #[derive(Default)]
    pub struct SimulationFlags: u32
    {
        /// If set, and BUFFER_EDGES is not set,
        /// then simplification will use a reusable set
        /// of buffers for each call.  Otherwise,
        /// these buffers will be allocated each time
        /// simplification happens.
        const USE_STATE = 1 << 0;
        /// If set, edge buffering will be used.
        /// If not set, then the standard "record
        /// and sort" method will be used.
        const BUFFER_EDGES = 1 << 1;
    }
}

#[derive(Copy, Clone, Debug)]
pub struct SimulationParams {
    pub popsize: u32,
    pub mutrate: f64,
    pub psurvival: f64,
    pub xovers: f64,
    pub genome_length: Position,
    pub buffer_edges: bool,
    pub simplification_interval: Option<i64>,
    pub seed: usize,
    pub nsteps: i64,
    pub flags: SimulationFlags,
    pub simplification_flags: SimplificationFlags,
}

impl Default for SimulationParams {
    fn default() -> Self {
        Self {
            popsize: 1000,
            mutrate: 0.0,
            psurvival: 0.0,
            xovers: 0.0,
            genome_length: 1000000.into(),
            buffer_edges: false,
            simplification_interval: None,
            seed: 0,
            nsteps: 0,
            flags: SimulationFlags::empty(),
            simplification_flags: SimplificationFlags::empty(),
        }
    }
}

fn mutate_tables(mutrate: f64, tables: &mut TableCollection, rng: &mut Rng) {
    match mutrate.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => (),
        Some(_) => return,
        None => panic!("bad mutation rate"),
    };
    let mut posmap = std::collections::HashMap::<PositionLLType, SiteId>::new();
    let mut derived_map = std::collections::HashMap::<PositionLLType, u8>::new();

    let num_edges = tables.edges().len();
    for i in 0..num_edges {
        let e = *tables.edge(EdgeId::from(i));
        let ptime = i64::from(tables.node(e.parent).time);
        let ctime = i64::from(tables.node(e.child).time);
        let blen = ctime - ptime;
        assert!((blen as i64) > 0, "{} {} {}", blen, ptime, ctime,);

        let pedge = ((PositionLLType::from(e.right) - PositionLLType::from(e.left)) as f64)
            / (PositionLLType::from(tables.genome_length()) as f64);

        let mutrate_edge = (mutrate * blen as f64) * pedge;
        let nmuts = rng.0.poisson(mutrate_edge);
        for _ in 0..nmuts {
            let t = ((rng.0.flat(ptime as f64, ctime as f64) as i64) + 1) as f64;
            let pos = rng.0.flat(
                PositionLLType::from(e.left) as f64,
                PositionLLType::from(e.right) as f64,
            ) as PositionLLType;

            match posmap.get(&pos) {
                Some(x) => {
                    // Get a new derived state for this site
                    let dstate = match derived_map.get(&pos) {
                        Some(y) => y + 1,
                        None => 1,
                    };
                    derived_map.insert(pos, dstate).unwrap();
                    tables
                        .add_mutation(e.child, None, *x, t, Some(vec![dstate]), true)
                        .unwrap();
                }
                None => {
                    let site_id = tables.add_site(pos, Some(vec![0])).unwrap();
                    tables
                        .add_mutation(e.child, None, site_id, t, Some(vec![1]), true)
                        .unwrap();

                    if posmap.insert(pos, site_id).is_some() {
                        panic!("hash failure");
                    }
                    if derived_map.insert(pos, 1).is_some() {
                        panic!("derived state hash failure");
                    }
                }
            }
        }
    }
    tables.sort_tables(TableSortingFlags::SKIP_EDGE_TABLE);
}

pub fn neutral_wf(
    params: SimulationParams,
) -> Result<(TableCollection, Vec<i32>), Box<dyn std::error::Error>> {
    // FIXME: gotta validate input params!

    let mut actual_simplification_interval: i64 = -1;

    match params.simplification_interval {
        None => (),
        Some(x) => actual_simplification_interval = validate_simplification_interval(x),
    }

    let mut rng = Rng::new(params.seed);

    let mut pop = PopulationState::new(params.genome_length);
    let mut samples: SamplesInfo = Default::default();

    // Record nodes for the first generation
    // Nodes will have birth time 0 in deme 0.
    for index in 0..params.popsize {
        let node0 = pop.tables.add_node(0_f64, 0).unwrap();
        let node1 = pop.tables.add_node(0_f64, 0).unwrap();
        pop.parents.push(Parent {
            index: index as usize,
            node0,
            node1,
        });
    }

    for i in 0..pop.tables.num_nodes() {
        samples.edge_buffer_founder_nodes.push(i.into());
    }

    let mut simplified = false;
    let mut state = SimplificationBuffers::new();

    let mut output = SimplificationOutput::new();

    let new_edge_handler = if params.flags.contains(SimulationFlags::BUFFER_EDGES) {
        buffer_edges
    } else {
        record_edges
    };

    let recrate = match params.xovers.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => Some(params.xovers),
        Some(std::cmp::Ordering::Equal) => None,
        Some(std::cmp::Ordering::Less) | None => panic!("invalid recombination rate"),
    };

    for birth_time in 1..(params.nsteps + 1) {
        deaths_and_parents(params.psurvival, &mut rng, &mut pop);
        generate_births(
            recrate,
            birth_time.into(),
            &mut rng,
            &mut pop,
            &new_edge_handler,
        );
        if actual_simplification_interval != -1 && birth_time % actual_simplification_interval == 0
        {
            simplify_and_remap_nodes(
                params.flags,
                params.simplification_flags,
                &mut samples,
                &mut state,
                &mut pop,
                &mut output,
            );
            simplified = true;
        } else {
            simplified = false;
        }
    }

    if !simplified && actual_simplification_interval != -1 {
        simplify_and_remap_nodes(
            params.flags,
            params.simplification_flags,
            &mut samples,
            &mut state,
            &mut pop,
            &mut output,
        );
    }

    let mut is_alive: Vec<i32> = vec![0; pop.tables.num_nodes()];

    for p in pop.parents {
        is_alive[usize::from(p.node0)] = 1;
        is_alive[usize::from(p.node1)] = 1;
    }

    mutate_tables(params.mutrate, &mut pop.tables, &mut rng);

    for s in pop.tables.sites() {
        match &s.ancestral_state {
            Some(x) => {
                assert_eq!(x.len(), 1);
                assert_eq!(x[0], 0);
            }
            None => panic!("ancestral_state is None"),
        };
    }
    for m in pop.tables.mutations() {
        match &m.derived_state {
            Some(x) => {
                assert_eq!(x.len(), 1);
                assert!(x[0] > 0);
            }
            None => panic!("derived_state is None"),
        };
    }

    Ok((pop.tables, is_alive))
}

// Below is code for simplifying in a separate thread via channels.
struct SimplificationRoundTripData {
    samples: SamplesInfo,
    edge_buffer: EdgeBuffer,
    state: SimplificationBuffers,
    output: SimplificationOutput,
}

impl SimplificationRoundTripData {
    fn new(
        samples: SamplesInfo,
        edge_buffer: EdgeBuffer,
        state: SimplificationBuffers,
        output: SimplificationOutput,
    ) -> Self {
        Self {
            samples,
            edge_buffer,
            state,
            output,
        }
    }
}

// Take ownership, simplify, return ownership.
fn simplify_from_edge_buffer_channel(
    flags: SimplificationFlags,
    inputs: SimplificationRoundTripData,
    tables: Arc<Mutex<TableCollection>>,
) -> Result<SimplificationRoundTripData, Box<dyn std::error::Error>> {
    let mut state = inputs.state;
    let mut edge_buffer = inputs.edge_buffer;
    let mut output = inputs.output;

    let mut t = tables.lock().unwrap();

    simplify_from_edge_buffer(
        &inputs.samples,
        flags,
        &mut state,
        &mut edge_buffer,
        &mut t,
        &mut output,
    )?;

    Ok(SimplificationRoundTripData::new(
        inputs.samples,
        edge_buffer,
        state,
        output,
    ))
}

fn generate_births_v2(
    breakpoint: BreakpointFunction,
    birth_time: Time,
    genome_length: Position,
    births: &[Birth],
    rng: &mut StdRng,
    parents: &mut [Parent],
    new_nodes: &mut NodeTable,
    new_edges: &mut EdgeTable,
    next_node_id: &mut TablesIdInteger,
) {
    for b in births {
        // Add the new nodes, but don't use them for recording yet
        new_nodes.push(Node {
            time: birth_time,
            deme: 0.into(),
            flags: 0,
        });
        new_nodes.push(Node {
            time: birth_time,
            deme: 0.into(),
            flags: 0,
        });

        let new_node_0 = NodeId::from(*next_node_id);
        let new_node_1 = NodeId::from(*next_node_id + 1);

        *next_node_id += 2;

        // replacing crossover_and_record_edges here...
        for (p, c) in [(b.parent0, new_node_0), (b.parent1, new_node_1)] {
            let mut pnodes = (p.node0, p.node1);
            mendel(&mut pnodes, rng);
            if let Some(exp) = breakpoint {
                let mut current_pos: PositionLLType = 0;
                loop {
                    let next_length = (rng.sample(exp) as PositionLLType) + 1;
                    if current_pos + next_length < genome_length {
                        new_edges.push(Edge {
                            parent: pnodes.0,
                            child: c,
                            left: current_pos.into(),
                            right: (current_pos + next_length).into(),
                        });
                        current_pos += next_length;
                        std::mem::swap(&mut pnodes.0, &mut pnodes.1);
                    } else {
                        new_edges.push(Edge {
                            parent: pnodes.0,
                            child: c,
                            left: current_pos.into(),
                            right: genome_length,
                        });
                        break;
                    }
                }
            } else {
                new_edges.push(Edge {
                    parent: pnodes.0,
                    child: c,
                    left: 0.into(),
                    right: genome_length,
                });
            }
        }
        parents[b.index].index = b.index;
        parents[b.index].node0 = new_node_0;
        parents[b.index].node1 = new_node_1;
    }
}

enum Simplifying {
    No((SamplesInfo, SimplificationBuffers, SimplificationOutput)),
    Yes(thread::JoinHandle<SimplificationRoundTripData>),
}

fn dispatch_simplification(
    birth_time: i64,
    pop: &mut PopulationState,
    new_nodes: &mut NodeTable,
    new_edges: &mut EdgeTable,
    first_child_node_after_last_simplification: TablesIdInteger,
    flags: SimplificationFlags,
    tables: Arc<Mutex<TableCollection>>,
    samples: SamplesInfo,
    state: SimplificationBuffers,
    output: SimplificationOutput,
) -> Simplifying {
    // If new nodes is empty, there's no work to be done
    // and we can return consumed stuff
    if new_nodes.is_empty() {
        Simplifying::No((samples, state, output))
    } else {
        //println!("Firing off some simplification at {}", birth_time);
        // Else, we have to do some moves of the big
        // data structures and return a JoinHandle
        let mut edge_buffer = EdgeBuffer::default();

        let mut samples = samples;

        {
            let mut t = tables.lock().unwrap();
            // Transfer our edges
            let num_nodes = t.nodes().len() as TablesIdInteger;
            for edge in new_edges.drain(0..) {
                let p = match edge.parent >= first_child_node_after_last_simplification {
                    false => TablesIdInteger::from(edge.parent),
                    true => {
                        //println!(
                        //    "parent mapping: {}, {} {}-> {}",
                        //    TablesIdInteger::from(edge.parent),
                        //    num_nodes,
                        //    first_child_node_after_last_simplification,
                        //    TablesIdInteger::from(edge.parent) + num_nodes
                        //        - first_child_node_after_last_simplification
                        //);
                        TablesIdInteger::from(edge.parent) + num_nodes
                            - first_child_node_after_last_simplification
                    }
                };

                // TODO: this can/should be dispatched to a thread.
                let c = match edge.child >= first_child_node_after_last_simplification {
                    false => TablesIdInteger::from(edge.child),
                    true => {
                        //println!(
                        //    "child mapping: {}, {} {}-> {}",
                        //    TablesIdInteger::from(edge.child),
                        //    num_nodes,
                        //    first_child_node_after_last_simplification,
                        //    TablesIdInteger::from(edge.child) + num_nodes
                        //        - first_child_node_after_last_simplification
                        //);
                        TablesIdInteger::from(edge.child) + num_nodes
                            - first_child_node_after_last_simplification
                    }
                };
                assert!(
                    (c as usize) < t.nodes().len() + new_nodes.len(),
                    "{} {} {}",
                    c,
                    t.nodes().len(),
                    new_nodes.len()
                );
                edge_buffer
                    .record_edge(p, c, edge.left, edge.right)
                    .unwrap();
            }
            assert!(new_edges.is_empty());
            //fill_samples(&pop.parents, &mut samples);
            samples.samples.clear();
            for p in pop.parents.iter_mut() {
                if p.node0 >= first_child_node_after_last_simplification {
                    p.node0 = (TablesIdInteger::from(p.node0) + num_nodes
                        - first_child_node_after_last_simplification)
                        .into();
                }
                if p.node1 >= first_child_node_after_last_simplification {
                    p.node1 = (TablesIdInteger::from(p.node1) + num_nodes
                        - first_child_node_after_last_simplification)
                        .into();
                }
                samples.samples.push(p.node0);
                samples.samples.push(p.node1);
            }
            assert_eq!(pop.parents.len() * 2, samples.samples.len());
            // transfer over our new nodes
            let mut node_table = t.dump_node_table();
            node_table.append(new_nodes);
            t.set_node_table(node_table);
        }

        samples.edge_buffer_founder_nodes.clear();
        for (i, p) in &mut pop.parents.iter_mut().enumerate() {
            p.node0 = NodeId::from(2 * i); // utput.idmap[usize::from(p.node0)];
            p.node1 = NodeId::from(2 * i + 1); // utput.idmap[usize::from(p.node0)];
            samples.edge_buffer_founder_nodes.push(p.node0);
            samples.edge_buffer_founder_nodes.push(p.node1);
        }

        Simplifying::Yes(thread::spawn(move || {
            // consume data
            let inputs = SimplificationRoundTripData::new(samples, edge_buffer, state, output);

            // send data to simplification
            let outputs = simplify_from_edge_buffer_channel(flags, inputs, tables).unwrap();

            outputs
        }))
    }
}

pub fn neutral_wf_simplify_separate_thread(
    params: SimulationParams,
) -> Result<(TableCollection, Vec<i32>), Box<dyn std::error::Error>> {
    // FIXME: gotta validate input params!
    // TODO: require a simplification interval > 0
    if !params.flags.contains(SimulationFlags::BUFFER_EDGES) {
        return Err(Box::new(SimulationError::ErrorMessage(
            "simulation using threads requires edge buffering".to_string(),
        )));
    }

    let mut actual_simplification_interval: i64 = -1;

    let breakpoint: BreakpointFunction = match params.xovers.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => Some(
            Exp::new(params.xovers / PositionLLType::from(params.genome_length) as f64).unwrap(),
        ),
        Some(_) => None,
        None => panic!("invalid xovers: {}", params.xovers),
    };

    match params.simplification_interval {
        None => (),
        Some(x) => actual_simplification_interval = validate_simplification_interval(x),
    }

    let mut rng = StdRng::seed_from_u64(params.seed);

    let mut pop = PopulationState::new(params.genome_length);
    let mut samples: SamplesInfo = Default::default();

    // Record nodes for the first generation
    // Nodes will have birth time 0 in deme 0.
    let mut next_node_id: TablesIdInteger = 0;
    let mut _tables = TableCollection::new(params.genome_length).unwrap();
    for index in 0..params.popsize {
        let node0 = _tables.add_node(0_f64, 0).unwrap();
        let node1 = _tables.add_node(0_f64, 0).unwrap();
        pop.parents.push(Parent {
            index: index as usize,
            node0,
            node1,
        });
        next_node_id += 2;
    }
    assert_eq!(next_node_id, _tables.nodes().len() as TablesIdInteger);
    let mut first_child_node_after_last_simplification = next_node_id;

    for i in 0.._tables.num_nodes() {
        samples.edge_buffer_founder_nodes.push(i.into());
    }

    let genome_length = _tables.genome_length();
    let tables = Arc::new(Mutex::new(_tables));
    let mut simplified = false;
    let mut state = SimplificationBuffers::new();
    let mut output = SimplificationOutput::new();

    let mut new_nodes = NodeTable::default();
    let mut new_edges = EdgeTable::default();

    let mut birth_time: i64 = 1;

    loop {
        // Step 1: check if there's work to simplify

        //println!("check if we simplify at {}", birth_time);

        let simplifying = dispatch_simplification(
            birth_time,
            &mut pop,
            &mut new_nodes,
            &mut new_edges,
            first_child_node_after_last_simplification,
            params.simplification_flags,
            tables.clone(),
            samples,
            state,
            output,
        );

        //    // Join our simplification thread handle, if
        //    // it exists

        //    fill_samples(&pop.parents, &mut samples);
        //    // transfer over our new nodes
        //    let mut node_table = pop.tables.dump_node_table();
        //    node_table.append(&mut new_nodes);
        //    pop.tables.set_node_table(node_table);

        //    // consume data
        //    let inputs = SimplificationRoundTripData::new(
        //        samples,
        //        pop.edge_buffer,
        //        pop.tables,
        //        state,
        //        output,
        //    );
        //    // send data to simplification
        //    let outputs = simplify_from_edge_buffer_channel(params.simplification_flags, inputs)?;
        //    // get our data back
        //    pop.edge_buffer = outputs.edge_buffer;
        //    pop.tables = outputs.tables;
        //    output = outputs.output;
        //    state = outputs.state;
        //    samples = outputs.samples;
        //    next_node_id = pop.tables.nodes().len() as TablesIdInteger;
        //    // remap parent nodes
        //    for p in &mut pop.parents {
        //        p.node0 = output.idmap[usize::from(p.node0)];
        //        p.node1 = output.idmap[usize::from(p.node1)];
        //        assert!(pop.tables.node(p.node0).flags & NodeFlags::IS_SAMPLE.bits() > 0);
        //    }

        //    // Track what (remapped) nodes are now alive.
        //    samples.edge_buffer_founder_nodes.clear();
        //    for p in &pop.parents {
        //        samples.edge_buffer_founder_nodes.push(p.node0);
        //        samples.edge_buffer_founder_nodes.push(p.node1);
        //    }
        //    simplified = true;
        //} else {
        //    simplified = false;
        //}

        //if birth_time > params.nsteps {
        //    break;
        //}

        // record new data while simplification is happening
        match simplifying {
            Simplifying::No(data) => {
                //println!("nope at time {}", i64::from(birth_time));
                simplified = false;
                samples = data.0;
                state = data.1;
                output = data.2;
                for _ in 1..(actual_simplification_interval + 1) {
                    deaths_and_parents(params.psurvival, &mut rng, &mut pop);
                    generate_births_v2(
                        breakpoint,
                        birth_time.into(),
                        genome_length,
                        &mut pop.births,
                        &mut rng,
                        &mut pop.parents,
                        &mut new_nodes,
                        &mut new_edges,
                        &mut next_node_id,
                    );

                    birth_time += 1;

                    // We may exit if the simplification interval
                    // and/or the nsteps is a "funny" value
                    if birth_time > params.nsteps {
                        //println!("breaking in ::No");
                        break;
                    }
                }
            }
            Simplifying::Yes(handle) => {
                //println!("wrapping up simplification at {}", i64::from(birth_time));
                let outputs = handle.join().unwrap();
                simplified = true;
                output = outputs.output;
                state = outputs.state;
                samples = outputs.samples;

                // FIXME: this is probalby our issue, causing "bad" things to
                // happen w/new node IDs?
                next_node_id = samples.samples.len() as TablesIdInteger;
                first_child_node_after_last_simplification = next_node_id;
                // println!(
                //     "{} {} {} {}",
                //     next_node_id,
                //     first_child_node_after_last_simplification,
                //     new_nodes.len(),
                //     new_edges.len()
                // );
                // remap parent nodes
                // FIXME NOTE TODO: fascinating--the idmap is coming back funky?
                // {
                //     let t = tables.lock().unwrap();
                //     for p in &mut pop.parents {
                //         p.node0 = output.idmap[usize::from(p.node0)];
                //         p.node1 = output.idmap[usize::from(p.node1)];
                //         assert!(t.node(p.node0).flags & NodeFlags::IS_SAMPLE.bits() > 0);
                //     }
                // }

                // // TODO: we can save a loop by merging the pushes into
                // // the previous loop
                // // Track what (remapped) nodes are now alive.
                // samples.edge_buffer_founder_nodes.clear();
                // for p in &pop.parents {
                //     samples.edge_buffer_founder_nodes.push(p.node0);
                //     samples.edge_buffer_founder_nodes.push(p.node1);
                // }
                if birth_time > params.nsteps {
                    //println!("breaking on ::Yes");
                    break;
                }
                //if !new_nodes.is_empty() {
            }
        }
    }

    // for birth_time in 1..(params.nsteps + 1) {
    //     deaths_and_parents(params.psurvival, &mut rng, &mut pop);
    //     generate_births_v2(
    //         breakpoint,
    //         birth_time.into(),
    //         pop.tables.genome_length(),
    //         &mut pop.births,
    //         &mut rng,
    //         &mut pop.parents,
    //         &mut new_nodes,
    //         &mut pop.edge_buffer,
    //         &mut next_node_id,
    //     );

    //     if actual_simplification_interval != -1 && birth_time % actual_simplification_interval == 0
    //     {
    //         fill_samples(&pop.parents, &mut samples);
    //         // transfer over our new nodes
    //         let mut node_table = pop.tables.dump_node_table();
    //         node_table.append(&mut new_nodes);
    //         pop.tables.set_node_table(node_table);

    //         // consume data
    //         let inputs = SimplificationRoundTripData::new(
    //             samples,
    //             pop.edge_buffer,
    //             pop.tables,
    //             state,
    //             output,
    //         );
    //         // send data to simplification
    //         let outputs = simplify_from_edge_buffer_channel(params.simplification_flags, inputs)?;
    //         // get our data back
    //         pop.edge_buffer = outputs.edge_buffer;
    //         pop.tables = outputs.tables;
    //         output = outputs.output;
    //         state = outputs.state;
    //         samples = outputs.samples;
    //         next_node_id = pop.tables.nodes().len() as TablesIdInteger;
    //         // remap parent nodes
    //         for p in &mut pop.parents {
    //             p.node0 = output.idmap[usize::from(p.node0)];
    //             p.node1 = output.idmap[usize::from(p.node1)];
    //             assert!(pop.tables.node(p.node0).flags & NodeFlags::IS_SAMPLE.bits() > 0);
    //         }

    //         // Track what (remapped) nodes are now alive.
    //         samples.edge_buffer_founder_nodes.clear();
    //         for p in &pop.parents {
    //             samples.edge_buffer_founder_nodes.push(p.node0);
    //             samples.edge_buffer_founder_nodes.push(p.node1);
    //         }
    //         simplified = true;
    //     } else {
    //         simplified = false;
    //     }
    // }

    // if !simplified && actual_simplification_interval != -1 {
    //     if !new_nodes.is_empty() {
    //         let mut node_table = tables.dump_node_table();
    //         node_table.append(&mut new_nodes);
    //         tables.set_node_table(node_table);
    //     }
    //     simplify_and_remap_nodes(
    //         params.flags,
    //         params.simplification_flags,
    //         &mut samples,
    //         &mut state,
    //         &mut pop,
    //         &mut output,
    //     );
    // }

    let mut return_tables = match Arc::try_unwrap(tables) {
        Ok(x) => match x.into_inner() {
            Ok(tables) => tables,
            Err(_) => panic!("poisoned mutex"),
        },
        Err(_) => panic!("multiple references to tables still in play!"),
    };

    let mut is_alive: Vec<i32> = vec![0; return_tables.num_nodes()];

    for p in pop.parents {
        is_alive[usize::from(p.node0)] = 1;
        is_alive[usize::from(p.node1)] = 1;
    }

    mutate_tables(params.mutrate, &mut return_tables, &mut rng);

    for s in return_tables.sites() {
        match &s.ancestral_state {
            Some(x) => {
                assert_eq!(x.len(), 1);
                assert_eq!(x[0], 0);
            }
            None => panic!("ancestral_state is None"),
        };
    }
    for m in return_tables.mutations() {
        match &m.derived_state {
            Some(x) => {
                assert_eq!(x.len(), 1);
                assert!(x[0] > 0);
            }
            None => panic!("derived_state is None"),
        };
    }

    Ok((return_tables, is_alive))
}
