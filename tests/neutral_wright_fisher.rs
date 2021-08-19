use bitflags::bitflags;
use forrustts::*;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Exp, Uniform};

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
    recorder: &impl Fn(NodeId, NodeId, (Position, Position), &mut TableCollection, &mut EdgeBuffer),
    rng: &mut StdRng,
    tables: &mut TableCollection,
    edge_buffer: &mut EdgeBuffer,
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
    recorder: &impl Fn(NodeId, NodeId, (Position, Position), &mut TableCollection, &mut EdgeBuffer),
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
    _: &mut TableCollection,
    buffer: &mut EdgeBuffer,
) {
    buffer
        .extend(parent.into(), Segment::new(span.0, span.1, child))
        .unwrap();
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

fn mutate_tables(mutrate: f64, tables: &mut TableCollection, rng: &mut StdRng) {
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
        let mutrate_edge = (mutrate * blen as f64)
            / (PositionLLType::from(e.right) - PositionLLType::from(e.left)) as f64;
        let exp = Exp::new(mutrate_edge).unwrap();
        let mut pos = PositionLLType::from(e.left) + (rng.sample(exp) as PositionLLType) + 1;
        let make_time = Uniform::new(ptime, ctime);
        let mut next_key: usize = 0;
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
                    derived_map.insert(pos, dstate).unwrap();
                    tables
                        .add_mutation(e.child, next_key, *x, t, Some(vec![dstate]), true)
                        .unwrap();
                }
                None => {
                    let site_id = tables.add_site(pos, Some(vec![0])).unwrap();
                    tables
                        .add_mutation(e.child, next_key, site_id, t, Some(vec![1]), true)
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
            next_key += 1;
        }
    }
    tables.sort_tables(TableSortingFlags::SKIP_EDGE_TABLE);
}

pub fn add_tskit_mutation_site_tables(
    tables: &TableCollection,
    g: Time,
    tskit_tables: &mut tskit::TableCollection,
) {
    for s in tables.sites() {
        tskit_tables
            .add_site(
                PositionLLType::from(s.position) as f64,
                match &s.ancestral_state {
                    Some(x) => Some(x),
                    None => panic!("expected ancestral_state"),
                },
            )
            .unwrap();
    }

    for m in tables.mutations() {
        let reverser = tskit_tools::simple_time_reverser(g);
        tskit_tables
            .add_mutation(
                tskit::tsk_id_t::from(m.site),
                tskit::tsk_id_t::from(m.node),
                tskit::TSK_NULL,
                reverser(m.time),
                Some(m.derived_state.as_ref().unwrap()),
            )
            .unwrap();
    }
}

pub fn neutral_wf(params: SimulationParams) -> Result<(TableCollection, Vec<i32>), ForrusttsError> {
    // FIXME: gotta validate input params!

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
