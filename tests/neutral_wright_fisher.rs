use bitflags::bitflags;
use forrustts::*;
use rand::SeedableRng;
use rand_distr::Poisson;

// Some of the material below seems like a candidate for a public API,
// but we need to decide here if this package should provide that.
// If so, then many of these types should not be here, as they have nothing
// to do with Wright-Fisher itself, and are instead more general.

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
            i64::from(self.left),
            i64::from(self.right)
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

fn deaths_and_parents<R: rand::Rng>(psurvival: f64, rng: &mut R, pop: &mut PopulationState) {
    let u = rand::distributions::Uniform::new(0.0, 1.0);
    let uparents = rand::distributions::Uniform::new(0_usize, pop.parents.len());
    pop.births.clear();
    for i in 0..pop.parents.len() {
        let x = rng.sample(u);
        match x.partial_cmp(&psurvival) {
            Some(std::cmp::Ordering::Greater) => {
                let random_index = rng.sample(uparents);
                let parent0 = pop.parents[random_index];
                let random_index = rng.sample(uparents);
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

fn mendel<R: rand::Rng>(pnodes: &mut (NodeId, NodeId), rng: &mut R) {
    let u = rand::distributions::Uniform::new(0.0, 1.0);
    let x: f64 = rng.sample(u);
    match x.partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Less) => {
            std::mem::swap(&mut pnodes.0, &mut pnodes.1);
        }
        Some(_) => (),
        None => panic!("Unexpected None"),
    }
}

fn crossover_breakpoints<R: rand::Rng>(
    recrate: Option<f64>,
    genome_length: Position,
    rng: &mut R,
) -> Option<Vec<Position>> {
    match recrate {
        Some(x) => {
            let poiss = Poisson::new(x).unwrap();
            let rpos = rand::distributions::Uniform::new(0, i64::from(genome_length));
            let n = rng.sample(poiss) as u32;
            match n > 0 {
                true => {
                    let mut rv = vec![];
                    for _ in 0..n {
                        rv.push(Position::from(rng.sample(rpos)));
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

fn crossover_and_record_edges<R: rand::Rng>(
    parent: Parent,
    child: NodeId,
    recrate: Option<f64>,
    recorder: &impl Fn(NodeId, NodeId, (Position, Position), &mut TableCollection, &mut EdgeBuffer),
    rng: &mut R,
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

fn generate_births<R: rand::Rng>(
    recrate: Option<f64>,
    birth_time: Time,
    rng: &mut R,
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
        p.node0 = output.idmap[usize::try_from(p.node0).unwrap()];
        p.node1 = output.idmap[usize::try_from(p.node1).unwrap()];
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
    pub seed: u64,
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

fn mutate_tables<R: rand::Rng>(mutrate: f64, tables: &mut TableCollection, rng: &mut R) {
    match mutrate.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => (),
        Some(_) => return,
        None => panic!("bad mutation rate"),
    };
    let mut posmap = std::collections::HashMap::<i64, SiteId>::new();
    let mut derived_map = std::collections::HashMap::<i64, u8>::new();

    let num_edges = tables.edges().len();
    for i in 0..num_edges {
        let e = *tables.edge(EdgeId::try_from(i).unwrap());
        let ptime = i64::from(tables.node(e.parent).time);
        let ctime = i64::from(tables.node(e.child).time);
        let blen = ctime - ptime;
        assert!((blen as i64) > 0, "{} {} {}", blen, ptime, ctime,);

        let pedge = ((i64::from(e.right) - i64::from(e.left)) as f64)
            / (i64::from(tables.genome_length()) as f64);

        let mutrate_edge = (mutrate * blen as f64) * pedge;
        let nmut_generator = Poisson::new(mutrate_edge).unwrap();
        let nmuts = rng.sample(nmut_generator) as u32;
        let pos_generator = rand::distributions::Uniform::new(e.left.raw(), e.right.raw());
        let time_generator = rand::distributions::Uniform::new(ptime, ctime);
        for _ in 0..nmuts {
            let t = rng.sample(time_generator) + 1;
            let pos = rng.sample(pos_generator);

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

    let mut rng = rand::rngs::StdRng::seed_from_u64(params.seed);
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
        samples
            .edge_buffer_founder_nodes
            .push(i.try_into().unwrap());
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
        is_alive[usize::try_from(p.node0).unwrap()] = 1;
        is_alive[usize::try_from(p.node1).unwrap()] = 1;
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