use crate::newtypes::EdgeId;
use crate::newtypes::NodeId;
use crate::newtypes::Position;
use crate::newtypes::PositionLLType;
use crate::newtypes::SiteId;
use crate::newtypes::Time;
use crate::traits::AncestryType;
use crate::ForrusttsError;
use bitflags::bitflags;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Exp, Uniform};
use std::convert::TryFrom;
use std::convert::TryInto;
use tskit::TableAccess;
// Some of the material below seems like a candidate for a public API,
// but we need to decide here if this package should provide that.
// If so, then many of these types should not be here, as they have nothing
// to do with Wright-Fisher itself, and are instead more general.

// Even though Position is an integer, we will use
// an exponential distribution to get the distance to
// the next crossover position.  The reason for this is
// that rand_distr::Geometric has really poor performance.
type BreakpointFunction = Option<Exp<f64>>;

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
    pub edge_buffer: crate::EdgeBuffer,
    pub tables: crate::TableCollection,
}

impl PopulationState {
    pub fn new(genome_length: Position) -> Self {
        PopulationState {
            parents: vec![],
            births: vec![],
            edge_buffer: crate::EdgeBuffer::new(),
            tables: crate::TableCollection::new(genome_length).unwrap(),
        }
    }
}

fn deaths_and_parents(psurvival: f64, rng: &mut StdRng, pop: &mut PopulationState) {
    pop.births.clear();
    let random_parents = Uniform::new(0_usize, pop.parents.len() as usize);
    for i in 0..pop.parents.len() {
        let x: f64 = rng.gen();
        match x.partial_cmp(&psurvival) {
            Some(std::cmp::Ordering::Greater) => {
                let parent0 = pop.parents[rng.sample(random_parents)];
                let parent1 = pop.parents[rng.sample(random_parents)];
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

fn mendel(pnodes: &mut (NodeId, NodeId), rng: &mut StdRng) {
    let x: f64 = rng.gen();
    match x.partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Less) => {
            std::mem::swap(&mut pnodes.0, &mut pnodes.1);
        }
        Some(_) => (),
        None => panic!("Unexpected None"),
    }
}

fn crossover_and_record_edges(
    parent: Parent,
    child: NodeId,
    breakpoint: BreakpointFunction,
    recorder: &impl Fn(
        NodeId,
        NodeId,
        (Position, Position),
        &mut crate::TableCollection,
        &mut crate::EdgeBuffer,
    ),
    rng: &mut StdRng,
    tables: &mut crate::TableCollection,
    edge_buffer: &mut crate::EdgeBuffer,
) {
    let mut pnodes = (parent.node0, parent.node1);
    mendel(&mut pnodes, rng);

    if let Some(exp) = breakpoint {
        let mut current_pos: PositionLLType = 0;
        loop {
            // TODO: gotta justify the next line...
            let next_length = (rng.sample(exp) as PositionLLType) + 1;
            assert!(next_length > 0);
            if current_pos + next_length < tables.genome_length() {
                recorder(
                    pnodes.0,
                    child,
                    (current_pos.into(), (current_pos + next_length).into()),
                    tables,
                    edge_buffer,
                );
                current_pos += next_length;
                std::mem::swap(&mut pnodes.0, &mut pnodes.1);
            } else {
                recorder(
                    pnodes.0,
                    child,
                    (current_pos.into(), tables.genome_length()),
                    tables,
                    edge_buffer,
                );

                break;
            }
        }
    } else {
        recorder(
            pnodes.0,
            child,
            (0.into(), tables.genome_length()),
            tables,
            edge_buffer,
        );
    }
}

fn generate_births(
    breakpoint: BreakpointFunction,
    birth_time: Time,
    rng: &mut StdRng,
    pop: &mut PopulationState,
    recorder: &impl Fn(
        NodeId,
        NodeId,
        (Position, Position),
        &mut crate::TableCollection,
        &mut crate::EdgeBuffer,
    ),
) {
    for b in &pop.births {
        // Record 2 new nodes
        let new_node_0: NodeId = pop.tables.add_node(birth_time, 0).unwrap();
        let new_node_1: NodeId = pop.tables.add_node(birth_time, 0).unwrap();

        crossover_and_record_edges(
            b.parent0,
            new_node_0,
            breakpoint,
            recorder,
            rng,
            &mut pop.tables,
            &mut pop.edge_buffer,
        );
        crossover_and_record_edges(
            b.parent1,
            new_node_1,
            breakpoint,
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
    _: &mut crate::TableCollection,
    buffer: &mut crate::EdgeBuffer,
) {
    buffer
        .extend(parent.value(), crate::Segment::new(span.0, span.1, child))
        .unwrap();
}

fn record_edges(
    parent: NodeId,
    child: NodeId,
    span: (Position, Position),
    tables: &mut crate::TableCollection,
    _: &mut crate::EdgeBuffer,
) {
    tables.add_edge(span.0, span.1, parent, child).unwrap();
}

fn fill_samples(parents: &[Parent], samples: &mut crate::SamplesInfo) {
    samples.samples.clear();
    for p in parents {
        samples.samples.push(p.node0);
        samples.samples.push(p.node1);
    }
}

fn sort_and_simplify(
    flags: SimulationFlags,
    simplification_flags: crate::SimplificationFlags,
    samples: &crate::SamplesInfo,
    state: &mut crate::SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut crate::SimplificationOutput,
) {
    if !flags.contains(SimulationFlags::BUFFER_EDGES) {
        pop.tables.sort_tables(crate::TableSortingFlags::empty());
        if flags.contains(SimulationFlags::USE_STATE) {
            crate::simplify_tables(
                samples,
                simplification_flags,
                state,
                &mut pop.tables,
                output,
            )
            .unwrap();
        } else {
            crate::simplify_tables_without_state(
                samples,
                simplification_flags,
                &mut pop.tables,
                output,
            )
            .unwrap();
        }
    } else {
        crate::simplify_from_edge_buffer(
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
    simplification_flags: crate::SimplificationFlags,
    samples: &mut crate::SamplesInfo,
    state: &mut crate::SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut crate::SimplificationOutput,
) {
    fill_samples(&pop.parents, samples);
    sort_and_simplify(flags, simplification_flags, samples, state, pop, output);

    for p in &mut pop.parents {
        p.node0 = output.idmap[usize::try_from(p.node0).unwrap()];
        p.node1 = output.idmap[usize::try_from(p.node1).unwrap()];
        assert!(pop.tables.node(p.node0).flags & crate::NodeFlags::IS_SAMPLE.bits() > 0);
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
        // If set, and BUFFER_EDGES is not set,
        // then simplification will use a reusable set
        // of buffers for each call.  Otherwise,
        // these buffers will be allocated each time
        // simplification happens.
        const USE_STATE = 1 << 0;
        // If set, edge buffering will be used.
        // If not set, then the standard "record
        // and sort" method will be used.
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
    pub simplification_flags: crate::SimplificationFlags,
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
            simplification_flags: crate::SimplificationFlags::empty(),
        }
    }
}

fn mutate_tables(mutrate: f64, tables: &mut crate::TableCollection, rng: &mut StdRng) -> Vec<Time> {
    match mutrate.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => (),
        Some(_) => return vec![],
        None => panic!("bad mutation rate"),
    };
    let mut posmap = std::collections::HashMap::<PositionLLType, SiteId>::new();
    let mut derived_map = std::collections::HashMap::<PositionLLType, u8>::new();

    let mut origin_times_init: Vec<(Time, SiteId)> = vec![];
    let num_edges = tables.edges().len();
    for i in 0..num_edges {
        let e = *tables.edge(EdgeId::try_from(i).unwrap());
        let ptime = tables.node(e.parent).time.value() as i64;
        let ctime = tables.node(e.child).time.value() as i64;
        let blen = ctime - ptime;
        assert!((blen as i64) > 0, "{} {} {}", blen, ptime, ctime,);
        let mutrate_edge = (mutrate * blen as f64) / (e.right.value() - e.left.value()) as f64;
        let exp = Exp::new(mutrate_edge).unwrap();
        let mut pos = e.left.value() + (rng.sample(exp) as PositionLLType) + 1;
        let make_time = Uniform::new(ptime, ctime);
        while pos < e.right {
            assert!(ctime > ptime);
            let t = rng.sample(make_time) + 1;
            assert!(t <= ctime);
            assert!(t > ptime);
            match posmap.get(&pos) {
                Some(x) => {
                    // Get a new derived state for this site
                    let dstate = match derived_map.get(&pos) {
                        Some(y) => y + 1,
                        None => 1,
                    };
                    origin_times_init.push((t.into(), *x));
                    derived_map.insert(pos, dstate).unwrap();
                    tables
                        .add_mutation(
                            e.child,
                            origin_times_init.len() - 1,
                            *x,
                            Some(vec![dstate]),
                            true,
                        )
                        .unwrap();
                }
                None => {
                    let site_id = tables.add_site(pos, Some(vec![0])).unwrap();
                    origin_times_init.push((t.into(), site_id));
                    tables
                        .add_mutation(
                            e.child,
                            origin_times_init.len() - 1,
                            site_id,
                            Some(vec![1]),
                            true,
                        )
                        .unwrap();

                    if posmap.insert(pos, site_id).is_some() {
                        panic!("hash failure");
                    }
                    if derived_map.insert(pos, 1).is_some() {
                        panic!("derived state hash failure");
                    }
                }
            }
            pos += (rng.sample(exp) as PositionLLType) + 1;
        }
    }
    assert_eq!(origin_times_init.len(), tables.mutations().len());
    assert!(posmap.len() == derived_map.len());
    origin_times_init.sort_by(|a, b| {
        let pa = tables.site(a.1.try_into().unwrap()).position;
        let pb = tables.site(b.1.try_into().unwrap()).position;
        pa.cmp(&pb)
    });
    tables.sort_tables(crate::TableSortingFlags::SKIP_EDGE_TABLE);
    let mut rv = vec![];
    for (i, _) in origin_times_init {
        rv.push(i);
    }
    rv
}

fn add_tskit_mutation_site_tables(
    tables: &crate::TableCollection,
    origin_times: &[Time],
    g: Time,
    tskit_tables: &mut tskit::TableCollection,
) {
    for s in tables.sites() {
        tskit_tables
            .add_site(
                s.position.value() as f64,
                match &s.ancestral_state {
                    Some(x) => Some(x),
                    None => panic!("expected ancestral_state"),
                },
            )
            .unwrap();
    }

    for (i, m) in tables.enumerate_mutations() {
        let reverser = crate::tskit_tools::simple_time_reverser(g);
        assert!(match reverser(origin_times[i])
            .partial_cmp(&tskit_tables.nodes().time(m.node.value()).unwrap())
        {
            Some(std::cmp::Ordering::Less) => false,
            Some(_) => true,
            None => panic!("bad ordering"),
        });
        tskit_tables
            .add_mutation(
                m.site.value() as tskit::tsk_id_t,
                m.node.value() as tskit::tsk_id_t,
                tskit::TSK_NULL,
                reverser(origin_times[i]),
                Some(m.derived_state.as_ref().unwrap()),
            )
            .unwrap();
    }
}

pub fn neutral_wf(
    params: SimulationParams,
) -> Result<(crate::TableCollection, Vec<i32>, Vec<Time>), ForrusttsError> {
    // FIXME: gotta validate input params!

    let mut actual_simplification_interval: i64 = -1;

    let breakpoint: BreakpointFunction = match params.xovers.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => {
            Some(Exp::new(params.xovers / params.genome_length.value() as f64).unwrap())
        }
        Some(_) => None,
        None => panic!("invalid xovers: {}", params.xovers),
    };

    match params.simplification_interval {
        None => (),
        Some(x) => actual_simplification_interval = validate_simplification_interval(x),
    }

    let mut rng = StdRng::seed_from_u64(params.seed);

    let mut pop = PopulationState::new(params.genome_length);
    let mut samples: crate::SamplesInfo = Default::default();

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
    let mut state = crate::SimplificationBuffers::new();

    let mut output = crate::SimplificationOutput::new();

    let new_edge_handler = if params.flags.contains(SimulationFlags::BUFFER_EDGES) {
        buffer_edges
    } else {
        record_edges
    };

    for birth_time in 1..(params.nsteps + 1) {
        deaths_and_parents(params.psurvival, &mut rng, &mut pop);
        generate_births(
            breakpoint,
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

    let mut is_alive: Vec<i32> = vec![0.try_into().unwrap(); pop.tables.num_nodes()];

    for p in pop.parents {
        is_alive[usize::try_from(p.node0).unwrap()] = 1;
        is_alive[usize::try_from(p.node1).unwrap()] = 1;
    }

    let origin_times = mutate_tables(params.mutrate, &mut pop.tables, &mut rng);

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

    Ok((pop.tables, is_alive, origin_times))
}

pub struct SimulatorIterator {
    rng: StdRng,
    params: SimulationParams,
    make_seed: Uniform<u64>,
    nreps: i32,
    rep: i32,
}

pub struct SimResults {
    pub tables: crate::TableCollection,
    pub tsk_tables: tskit::TableCollection,
    pub is_sample: Vec<i32>,
}

impl SimulatorIterator {
    fn new(params: SimulationParams, nreps: i32) -> Self {
        Self {
            rng: StdRng::seed_from_u64(params.seed),
            params,
            make_seed: Uniform::new(0_u64, u32::MAX as u64),
            nreps,
            rep: 0,
        }
    }
}

impl Iterator for SimulatorIterator {
    type Item = SimResults;

    fn next(&mut self) -> Option<Self::Item> {
        if self.rep < self.nreps {
            let seed = self.rng.sample(self.make_seed);
            let mut params = self.params;
            params.seed = seed;
            let (mut tables, is_sample, origin_times) = neutral_wf(params).unwrap();

            tables.sort_tables(crate::TableSortingFlags::empty());

            let mut tsk_tables = crate::tskit_tools::convert_to_tskit_minimal(
                &tables,
                &is_sample,
                crate::tskit_tools::simple_time_reverser(params.nsteps.into()),
                // Do not index tables here!
                // Things are unsorted!
                false,
            );
            add_tskit_mutation_site_tables(
                &tables,
                &origin_times,
                params.nsteps.into(),
                &mut tsk_tables,
            );
            tsk_tables
                .full_sort(tskit::TableSortOptions::empty())
                .unwrap();
            self.rep += 1;
            Some(SimResults {
                tables,
                tsk_tables,
                is_sample,
            })
        } else {
            None
        }
    }
}

pub struct Simulator {
    params: SimulationParams,
    nreps: i32,
}

impl Simulator {
    pub fn new(params: SimulationParams, nreps: i32) -> Self {
        Self { params, nreps }
    }

    pub fn iter(&self) -> SimulatorIterator {
        SimulatorIterator::new(self.params, self.nreps)
    }
}

pub fn make_samples(l: &[i32]) -> crate::SamplesInfo {
    let mut rv = crate::SamplesInfo::default();
    for (i, j) in l.iter().enumerate() {
        if *j == 1 {
            rv.samples.push(i.try_into().unwrap());
        }
    }
    rv
}
