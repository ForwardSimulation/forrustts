use clap::{value_t, value_t_or_exit, App, Arg};
use forrustts::EdgeId;
use forrustts::ForrusttsError;
use forrustts::NodeId;
use forrustts::Position;
use forrustts::SiteId;
use forrustts::TableTypeIntoRaw;
use forrustts::Time;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Exp, Uniform};
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
    pub edge_buffer: forrustts::EdgeBuffer,
    pub tables: forrustts::TableCollection,
}

impl PopulationState {
    pub fn new(genome_length: Position) -> Self {
        PopulationState {
            parents: vec![],
            births: vec![],
            edge_buffer: forrustts::EdgeBuffer::new(),
            tables: forrustts::TableCollection::new(genome_length).unwrap(),
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
    rng: &mut StdRng,
    tables: &mut forrustts::TableCollection,
    edge_buffer: &mut forrustts::EdgeBuffer,
) {
    let mut pnodes = (parent.node0, parent.node1);
    mendel(&mut pnodes, rng);

    if let Some(exp) = breakpoint {
        let mut current_pos: forrustts::PositionLLType = 0;
        loop {
            // TODO: gotta justify the next line...
            let next_length = (rng.sample(exp) as forrustts::PositionLLType) + 1;
            assert!(next_length > 0);
            if current_pos + next_length < tables.genome_length() {
                buffer_edges(
                    pnodes.0,
                    child,
                    (current_pos, current_pos + next_length),
                    tables,
                    edge_buffer,
                );
                current_pos += next_length;
                std::mem::swap(&mut pnodes.0, &mut pnodes.1);
            } else {
                buffer_edges(
                    pnodes.0,
                    child,
                    (current_pos, tables.genome_length()),
                    tables,
                    edge_buffer,
                );

                break;
            }
        }
    } else {
        buffer_edges(
            pnodes.0,
            child,
            (0, tables.genome_length()),
            tables,
            edge_buffer,
        );
    }
}

fn generate_births(
    breakpoint: BreakpointFunction,
    birth_time: i64,
    rng: &mut StdRng,
    pop: &mut PopulationState,
) {
    for b in &pop.births {
        // Record 2 new nodes
        let new_node_0: NodeId = pop.tables.add_node(birth_time, 0).unwrap();
        let new_node_1: NodeId = pop.tables.add_node(birth_time, 0).unwrap();

        crossover_and_record_edges(
            b.parent0,
            new_node_0,
            breakpoint,
            rng,
            &mut pop.tables,
            &mut pop.edge_buffer,
        );
        crossover_and_record_edges(
            b.parent1,
            new_node_1,
            breakpoint,
            rng,
            &mut pop.tables,
            &mut pop.edge_buffer,
        );

        pop.parents[b.index].index = b.index;
        pop.parents[b.index].node0 = new_node_0;
        pop.parents[b.index].node1 = new_node_1;
    }
}

fn buffer_edges<P: Into<forrustts::Position>, P2: Into<forrustts::Position>>(
    parent: NodeId,
    child: NodeId,
    span: (P, P2),
    _: &mut forrustts::TableCollection,
    buffer: &mut forrustts::EdgeBuffer,
) {
    buffer
        .extend(
            parent.into(),
            forrustts::Segment::new(span.0.into(), span.1.into(), child),
        )
        .unwrap();
}

fn fill_samples(parents: &[Parent], samples: &mut forrustts::SamplesInfo) {
    samples.samples.clear();
    for p in parents {
        samples.samples.push(p.node0);
        samples.samples.push(p.node1);
    }
}

fn simplify(
    simplification_flags: forrustts::SimplificationFlags,
    samples: &forrustts::SamplesInfo,
    state: &mut forrustts::SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut forrustts::SimplificationOutput,
) {
    forrustts::simplify_from_edge_buffer(
        samples,
        simplification_flags,
        state,
        &mut pop.edge_buffer,
        &mut pop.tables,
        output,
    )
    .unwrap();
}

fn simplify_and_remap_nodes(
    simplification_flags: forrustts::SimplificationFlags,
    samples: &mut forrustts::SamplesInfo,
    state: &mut forrustts::SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut forrustts::SimplificationOutput,
) {
    fill_samples(&pop.parents, samples);
    simplify(simplification_flags, samples, state, pop, output);

    for p in &mut pop.parents {
        p.node0 = output.idmap[usize::from(p.node0)];
        p.node1 = output.idmap[usize::from(p.node1)];
        assert!(pop.tables.node(p.node0).flags & forrustts::NodeFlags::IS_SAMPLE.bits() > 0);
    }

    samples.edge_buffer_founder_nodes.clear();
    for p in &pop.parents {
        samples.edge_buffer_founder_nodes.push(p.node0);
        samples.edge_buffer_founder_nodes.push(p.node1);
    }
}

fn validate_simplification_interval(x: i64) -> i64 {
    if x < 1 {
        panic!("simplification_interval must be None or >= 1");
    }
    x
}

pub struct PopulationParams {
    pub size: u32,
    pub genome_length: Position,
    pub breakpoint: BreakpointFunction,
    pub psurvival: f64,
    pub mutrate: f64,
}

impl PopulationParams {
    pub fn new(
        size: u32,
        genome_length: Position,
        xovers: f64,
        psurvival: f64,
        mutrate: f64,
    ) -> Self {
        PopulationParams {
            size,
            genome_length,
            breakpoint: match xovers.partial_cmp(&0.0) {
                Some(std::cmp::Ordering::Greater) => {
                    Some(Exp::new(xovers / genome_length.into_raw() as f64).unwrap())
                }
                Some(_) => None,
                None => None,
            },
            psurvival,
            mutrate,
        }
    }
}

pub struct SimulationParams {
    pub simplification_interval: Option<i64>,
    pub seed: u64,
    pub nsteps: i64,
    pub simplification_flags: forrustts::SimplificationFlags,
}

impl SimulationParams {
    pub fn new(simplification_interval: Option<i64>, seed: u64, nsteps: i64) -> Self {
        SimulationParams {
            simplification_interval,
            seed,
            nsteps,
            simplification_flags: forrustts::SimplificationFlags::empty(),
        }
    }
}

fn mutate_tables(
    mutrate: f64,
    tables: &mut forrustts::TableCollection,
    rng: &mut StdRng,
) -> Vec<Time> {
    match mutrate.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => (),
        Some(_) => return vec![],
        None => panic!("bad mutation rate"),
    };
    let mut posmap = std::collections::HashMap::<forrustts::PositionLLType, SiteId>::new();
    let mut derived_map = std::collections::HashMap::<forrustts::PositionLLType, u8>::new();

    let mut origin_times_init: Vec<(Time, SiteId)> = vec![];
    let num_edges = tables.edges().len();
    for i in 0..num_edges {
        let e = *tables.edge(EdgeId::from(i));
        let ptime = i64::from(tables.node(e.parent).time);
        let ctime = i64::from(tables.node(e.child).time);
        let blen = ctime - ptime;
        assert!((blen as i64) > 0, "{} {} {}", blen, ptime, ctime,);
        let mutrate_edge =
            (mutrate * blen as f64) / (e.right.into_raw() - e.left.into_raw()) as f64;
        let exp = Exp::new(mutrate_edge).unwrap();
        let mut pos = e.left.into_raw() + (rng.sample(exp) as forrustts::PositionLLType) + 1;
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
                            t,
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
                            t,
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
            pos += (rng.sample(exp) as forrustts::PositionLLType) + 1;
        }
    }
    assert_eq!(origin_times_init.len(), tables.mutations().len());
    assert!(posmap.len() == derived_map.len());
    origin_times_init.sort_by(|a, b| {
        let pa = tables.site(a.1).position;
        let pb = tables.site(b.1).position;
        pa.cmp(&pb)
    });
    tables.sort_tables(forrustts::TableSortingFlags::SKIP_EDGE_TABLE);
    let mut rv = vec![];
    for (i, _) in origin_times_init {
        rv.push(i);
    }
    rv
}

fn add_tskit_mutation_site_tables(
    tables: &forrustts::TableCollection,
    origin_times: &[Time],
    g: Time,
    tskit_tables: &mut tskit::TableCollection,
) {
    for s in tables.sites() {
        tskit_tables
            .add_site(
                s.position.into_raw() as f64,
                match &s.ancestral_state {
                    Some(x) => Some(x),
                    None => panic!("expected ancestral_state"),
                },
            )
            .unwrap();
    }

    for (i, m) in tables.enumerate_mutations() {
        let reverser = forrustts::tskit_tools::simple_time_reverser(g);
        assert!(match reverser(origin_times[i])
            .partial_cmp(&tskit_tables.nodes().time(m.node.into()).unwrap())
        {
            Some(std::cmp::Ordering::Less) => false,
            Some(_) => true,
            None => panic!("bad ordering"),
        });
        tskit_tables
            .add_mutation(
                i32::from(m.site) as tskit::tsk_id_t,
                i32::from(m.node) as tskit::tsk_id_t,
                tskit::TSK_NULL,
                reverser(origin_times[i]),
                Some(m.derived_state.as_ref().unwrap()),
            )
            .unwrap();
    }
}

pub fn neutral_wf(
    pop_params: PopulationParams,
    params: SimulationParams,
) -> Result<(forrustts::TableCollection, Vec<i32>, Vec<Time>), ForrusttsError> {
    // FIXME: gotta validate input params!

    let mut actual_simplification_interval: i64 = -1;

    match params.simplification_interval {
        None => (),
        Some(x) => actual_simplification_interval = validate_simplification_interval(x),
    }

    let mut rng = StdRng::seed_from_u64(params.seed);

    let mut pop = PopulationState::new(pop_params.genome_length);
    let mut samples: forrustts::SamplesInfo = Default::default();

    // Record nodes for the first generation
    // Nodes will have birth time 0 in deme 0.
    for index in 0..pop_params.size {
        let node0 = pop.tables.add_node(0., 0).unwrap();
        let node1 = pop.tables.add_node(0., 0).unwrap();
        pop.parents.push(Parent {
            index: index as usize,
            node0,
            node1,
        });
    }

    for i in 0..pop.tables.num_nodes() {
        samples.edge_buffer_founder_nodes.push(NodeId::from(i));
    }

    let mut simplified = false;
    let mut state = forrustts::SimplificationBuffers::new();

    let mut output = forrustts::SimplificationOutput::new();

    for birth_time in 1..(params.nsteps + 1) {
        deaths_and_parents(pop_params.psurvival, &mut rng, &mut pop);
        generate_births(pop_params.breakpoint, birth_time, &mut rng, &mut pop);
        if actual_simplification_interval != -1 && birth_time % actual_simplification_interval == 0
        {
            simplify_and_remap_nodes(
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
            params.simplification_flags,
            &mut samples,
            &mut state,
            &mut pop,
            &mut output,
        );
    } else if actual_simplification_interval == -1 {
        // We've done a big sim, but never simplified.
        // This, our edges are all "trapped" in the edge buffer.
        // We need to copy them over, which we will do so
        // that the output is sorted by birth time.
        let mut edges = pop.tables.dump_edges();
        assert!(edges.is_empty());
        for (i, h) in pop.edge_buffer.head_iter().rev().enumerate() {
            if *h != forrustts::nested_forward_list::NULL_INDEX {
                let head = pop.edge_buffer.len() - i - 1;
                assert!(
                    (head as usize) < pop.tables.num_nodes(),
                    "{} {}",
                    *h,
                    pop.tables.num_nodes()
                );
                for seg in pop
                    .edge_buffer
                    .values_iter(head as forrustts::nested_forward_list::IndexType)
                {
                    assert!(usize::from(seg.node) < pop.tables.num_nodes());
                    edges.push(forrustts::Edge {
                        left: seg.left,
                        right: seg.right,
                        parent: head.into(),
                        child: seg.node,
                    });
                }
            }
        }
        pop.tables.set_edge_table(edges);
    }

    let mut is_alive: Vec<i32> = vec![0; pop.tables.num_nodes()];

    for p in pop.parents {
        is_alive[usize::from(p.node0)] = 1;
        is_alive[usize::from(p.node1)] = 1;
    }

    let origin_times = mutate_tables(pop_params.mutrate, &mut pop.tables, &mut rng);

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

fn main() {
    let matches = App::new("forward_simulation")
        .arg(
            Arg::with_name("popsize")
                .short("N")
                .long("popsize")
                .help("Diploid population size")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("nsteps")
                .short("n")
                .long("nsteps")
                .help("number of birth steps to simulate")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("xovers")
                .short("x")
                .long("xovers")
                .help("Mean number of crossovers per meiosis.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mutrate")
                .short("m")
                .long("mutrate")
                .help("Mean number of mutations per new gamete.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("simplification_interval")
                .short("s")
                .long("simplify")
                .help("number of generations between simplifications")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("outfile")
                .short("o")
                .long("outfile")
                .help("Name of output file. The format is a tskit \"trees\" file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("seed")
                .short("S")
                .long("seed")
                .help("Random number seed")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("buffer_edges")
                .short("B")
                .long("buffer_edges")
                .help("Use edge buffering instead of sorting")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("psurvival")
                .short("P")
                .long("psurvival")
                .help("Survival probability")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("validate_tables")
                .short("v")
                .long("validate_tables")
                .help("Validate all tables prior to simplification")
                .takes_value(false),
        )
        .get_matches();

    // TODO: default params

    let popsize = value_t_or_exit!(matches.value_of("popsize"), u32);
    let g = value_t_or_exit!(matches.value_of("nsteps"), i64);
    let xovers = value_t_or_exit!(matches.value_of("xovers"), f64);
    let mutrate = value_t_or_exit!(matches.value_of("mutrate"), f64);
    let simplify_input = value_t!(matches.value_of("simplification_interval"), i64).unwrap_or(-1);
    let psurvival = value_t!(matches.value_of("psurvival"), f64).unwrap_or(0.0);
    let seed = value_t_or_exit!(matches.value_of("seed"), u64);
    let outfile = value_t_or_exit!(matches.value_of("outfile"), String);
    let validate_tables = matches.is_present("validate_tables");

    // TODO: parameter validation..

    let mut simplify: Option<i64> = None;

    if simplify_input > 0 {
        simplify = Some(simplify_input);
    }

    let mut simplification_flags = forrustts::SimplificationFlags::empty();

    if validate_tables {
        simplification_flags |= forrustts::SimplificationFlags::VALIDATE_ALL;
    }

    let (mut tables, is_sample, origin_times) = neutral_wf(
        PopulationParams::new(popsize, 10000000.into(), xovers, psurvival, mutrate),
        SimulationParams {
            simplification_interval: simplify,
            seed,
            nsteps: g,
            simplification_flags,
        },
    )
    .unwrap();

    let mut tskit_tables = forrustts::tskit_tools::convert_to_tskit_and_drain_minimal(
        &is_sample,
        forrustts::tskit_tools::simple_time_reverser(g.into()),
        simplify.is_some(),
        &mut tables,
    );

    match simplify.is_some() {
        true => (),
        false => {
            let _ = tskit_tables.build_index().unwrap();
        }
    };

    add_tskit_mutation_site_tables(&tables, &origin_times, g.into(), &mut tskit_tables);

    tskit_tables
        .dump(&outfile, tskit::TableOutputOptions::default())
        .unwrap();
}
