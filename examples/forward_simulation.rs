use bitflags::bitflags;
use clap::{value_t, value_t_or_exit, App, Arg};
use forrustts::ForrusttsError;
use forrustts::IdType;
use forrustts::Position;
use forrustts::Time;
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
    node0: IdType,
    node1: IdType,
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

fn mendel(pnodes: &mut (tskit::tsk_id_t, tskit::tsk_id_t), rng: &mut StdRng) {
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
    child: IdType,
    breakpoint: BreakpointFunction,
    recorder: &impl Fn(
        IdType,
        IdType,
        (Position, Position),
        &mut forrustts::TableCollection,
        &mut forrustts::EdgeBuffer,
    ),
    rng: &mut StdRng,
    tables: &mut forrustts::TableCollection,
    edge_buffer: &mut forrustts::EdgeBuffer,
) {
    let mut pnodes = (parent.node0, parent.node1);
    mendel(&mut pnodes, rng);
    let mut p0 = parent.node0;
    let mut p1 = parent.node1;

    if let Some(exp) = breakpoint {
        let mut current_pos: Position = 0;
        loop {
            // TODO: gotta justify the next line...
            let next_length = (rng.sample(exp) as Position) + 1;
            assert!(next_length > 0);
            if current_pos + next_length < tables.genome_length() {
                recorder(
                    p0,
                    child,
                    (current_pos, current_pos + next_length),
                    tables,
                    edge_buffer,
                );
                current_pos += next_length;
                std::mem::swap(&mut p0, &mut p1);
            } else {
                recorder(
                    p0,
                    child,
                    (current_pos, tables.genome_length()),
                    tables,
                    edge_buffer,
                );

                break;
            }
        }
    } else {
        recorder(p0, child, (0, tables.genome_length()), tables, edge_buffer);
    }
}

fn generate_births(
    breakpoint: BreakpointFunction,
    birth_time: Time,
    rng: &mut StdRng,
    pop: &mut PopulationState,
    recorder: &impl Fn(
        IdType,
        IdType,
        (Position, Position),
        &mut forrustts::TableCollection,
        &mut forrustts::EdgeBuffer,
    ),
) {
    for b in &pop.births {
        // Record 2 new nodes
        let new_node_0: IdType = pop.tables.add_node(birth_time, 0).unwrap();
        let new_node_1: IdType = pop.tables.add_node(birth_time, 0).unwrap();

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
    parent: IdType,
    child: IdType,
    span: (Position, Position),
    _: &mut forrustts::TableCollection,
    buffer: &mut forrustts::EdgeBuffer,
) {
    buffer
        .extend(parent, forrustts::Segment::new(span.0, span.1, child))
        .unwrap();
}

fn record_edges(
    parent: IdType,
    child: IdType,
    span: (Position, Position),
    tables: &mut forrustts::TableCollection,
    _: &mut forrustts::EdgeBuffer,
) {
    tables.add_edge(span.0, span.1, parent, child).unwrap();
}

fn fill_samples(parents: &[Parent], samples: &mut forrustts::SamplesInfo) {
    samples.samples.clear();
    for p in parents {
        samples.samples.push(p.node0);
        samples.samples.push(p.node1);
    }
}

fn sort_and_simplify(
    flags: SimulationFlags,
    simplification_flags: forrustts::SimplificationFlags,
    samples: &forrustts::SamplesInfo,
    state: &mut forrustts::SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut forrustts::SimplificationOutput,
) {
    if !flags.contains(SimulationFlags::BUFFER_EDGES) {
        pop.tables
            .sort_tables(forrustts::TableSortingFlags::empty());
        if flags.contains(SimulationFlags::USE_STATE) {
            forrustts::simplify_tables(
                samples,
                simplification_flags,
                state,
                &mut pop.tables,
                output,
            )
            .unwrap();
        } else {
            forrustts::simplify_tables_without_state(
                samples,
                simplification_flags,
                &mut pop.tables,
                output,
            )
            .unwrap();
        }
    } else {
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
}

fn simplify_and_remap_nodes(
    flags: SimulationFlags,
    simplification_flags: forrustts::SimplificationFlags,
    samples: &mut forrustts::SamplesInfo,
    state: &mut forrustts::SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut forrustts::SimplificationOutput,
) {
    fill_samples(&pop.parents, samples);
    sort_and_simplify(flags, simplification_flags, samples, state, pop, output);

    for p in &mut pop.parents {
        p.node0 = output.idmap[p.node0 as usize];
        p.node1 = output.idmap[p.node1 as usize];
    }

    if flags.contains(SimulationFlags::BUFFER_EDGES) {
        samples.edge_buffer_founder_nodes.clear();
        for p in &pop.parents {
            samples.edge_buffer_founder_nodes.push(p.node0);
            samples.edge_buffer_founder_nodes.push(p.node1);
        }
    }
}

fn validate_simplification_interval(x: Time) -> Time {
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
}

impl PopulationParams {
    pub fn new(size: u32, genome_length: Position, xovers: f64, psurvival: f64) -> Self {
        PopulationParams {
            size,
            genome_length,
            breakpoint: match xovers.partial_cmp(&0.0) {
                Some(std::cmp::Ordering::Greater) => {
                    Some(Exp::new(xovers / genome_length as f64).unwrap())
                }
                Some(_) => None,
                None => None,
            },
            psurvival,
        }
    }
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

pub struct SimulationParams {
    pub simplification_interval: Option<Time>,
    pub seed: u64,
    pub nsteps: Time,
    pub flags: SimulationFlags,
    pub simplification_flags: forrustts::SimplificationFlags,
}

impl SimulationParams {
    pub fn new(
        simplification_interval: Option<Time>,
        seed: u64,
        nsteps: Time,
        flags: SimulationFlags,
    ) -> Self {
        SimulationParams {
            simplification_interval,
            seed,
            nsteps,
            flags,
            simplification_flags: forrustts::SimplificationFlags::empty(),
        }
    }
}

pub fn neutral_wf(
    pop_params: PopulationParams,
    params: SimulationParams,
) -> Result<(forrustts::TableCollection, Vec<i32>), ForrusttsError> {
    // FIXME: gotta validate input params!

    let mut actual_simplification_interval: Time = -1;

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
        let node0 = pop.tables.add_node(0, 0).unwrap();
        let node1 = pop.tables.add_node(0, 0).unwrap();
        pop.parents.push(Parent {
            index: index as usize,
            node0,
            node1,
        });
    }

    for i in 0..pop.tables.num_nodes() {
        samples.edge_buffer_founder_nodes.push(i as IdType);
    }

    let mut simplified = false;
    let mut state = forrustts::SimplificationBuffers::new();

    let mut output = forrustts::SimplificationOutput::new();

    let new_edge_handler = if params.flags.contains(SimulationFlags::BUFFER_EDGES) {
        buffer_edges
    } else {
        record_edges
    };

    for birth_time in 1..(params.nsteps + 1) {
        deaths_and_parents(pop_params.psurvival, &mut rng, &mut pop);
        generate_births(
            pop_params.breakpoint,
            birth_time,
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
        is_alive[p.node0 as usize] = 1;
        is_alive[p.node1 as usize] = 1;
    }

    Ok((pop.tables, is_alive))
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

    let flags = if matches.is_present("buffer_edges") {
        SimulationFlags::BUFFER_EDGES
    } else {
        SimulationFlags::USE_STATE
    };

    let mut simplification_flags = forrustts::SimplificationFlags::empty();

    if validate_tables {
        simplification_flags |= forrustts::SimplificationFlags::VALIDATE_ALL;
    }

    let (mut tables, is_sample) = neutral_wf(
        PopulationParams::new(popsize, 10000000, xovers, psurvival),
        SimulationParams {
            simplification_interval: simplify,
            seed,
            nsteps: g,
            flags,
            simplification_flags,
        },
    )
    .unwrap();

    let mut tskit_tables = forrustts::tskit_tools::convert_to_tskit_and_drain_minimal(
        &is_sample,
        forrustts::tskit_tools::simple_time_reverser(g),
        simplify.is_some(),
        &mut tables,
    );

    match simplify.is_some() {
        true => (),
        false => {
            let _ = tskit_tables.build_index().unwrap();
        }
    };

    tskit_tables
        .dump(&outfile, tskit::TableOutputOptions::default())
        .unwrap();
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_flags() {
        let flags = SimulationFlags::empty();
        assert!(!flags.contains(SimulationFlags::USE_STATE));
        assert!(!flags.contains(SimulationFlags::BUFFER_EDGES));
        let flags = SimulationFlags::USE_STATE;
        assert!(flags.contains(SimulationFlags::USE_STATE));
        assert!(!flags.contains(SimulationFlags::BUFFER_EDGES));
        let flags = SimulationFlags::BUFFER_EDGES;
        assert!(!flags.contains(SimulationFlags::USE_STATE));
        assert!(flags.contains(SimulationFlags::BUFFER_EDGES));
    }

    fn run_sim(use_state: bool) -> (forrustts::TableCollection, Vec<i32>) {
        let flags = if use_state {
            SimulationFlags::USE_STATE
        } else {
            SimulationFlags::empty()
        };
        neutral_wf(
            PopulationParams::new(1000, 100000, 5e-3, 0.0),
            SimulationParams::new(Some(100), 666, 2000, flags),
        )
        .unwrap()
    }

    #[test]
    fn compare_state_to_no_state() {
        let (tables, _) = run_sim(false);
        let (tables_state, _) = run_sim(true);

        assert_eq!(tables.num_nodes(), tables_state.num_nodes());
        assert_eq!(tables.num_edges(), tables_state.num_edges());

        for (i, j) in tables.nodes().iter().zip(tables_state.nodes()) {
            assert_eq!(i.time, j.time);
            assert_eq!(i.deme, j.deme);
        }
        for (i, j) in tables.edges().iter().zip(tables_state.edges()) {
            assert_eq!(i.left, j.left);
            assert_eq!(i.right, j.right);
            assert_eq!(i.parent, j.parent);
            assert_eq!(i.child, j.child);
        }
    }
}

#[cfg(test)]
mod test_simplify_tables {
    use super::*;
    use forrustts::simplify_tables_without_state;
    use forrustts::ForrusttsError;
    use forrustts::SamplesInfo;
    use forrustts::SimplificationFlags;
    use forrustts::SimplificationOutput;
    use forrustts::TableCollection;
    use forrustts::{IdType, Position, Time};
    use tskit::TableAccess;

    fn simulate_data(
        num_generations: Time,
        genome_length: Position,
        psurvival: f64,
        seed: u64,
        // None here means "never simplify".
        simplification_interval: Option<Time>,
        flags: SimulationFlags,
    ) -> Result<(TableCollection, Vec<i32>), ForrusttsError> {
        neutral_wf(
            PopulationParams::new(250, genome_length, 5e-3, psurvival),
            SimulationParams {
                simplification_interval,
                seed,
                nsteps: num_generations,
                flags,
                simplification_flags: SimplificationFlags::VALIDATE_ALL,
            },
        )
    }

    #[test]
    fn test_kc_distance_to_tskit() {
        // 1. Simulate in rust
        // 2. Copy to tskit and reverse time
        // 3. sort and simplify via tskit
        // 4. sort and simplify via rust
        // 5. Do we get the same stuff out?

        let num_generations = 5000;
        let genome_length = 1000000;

        let (mut tables, mut is_sample) = simulate_data(
            num_generations,
            genome_length,
            0.0,
            42,
            None,
            SimulationFlags::empty(),
        )
        .unwrap();

        let mut tsk_tables = forrustts::tskit_tools::convert_to_tskit_minimal(
            &tables,
            &is_sample,
            forrustts::tskit_tools::simple_time_reverser(num_generations),
            // Do not index tables here!
            // Things are unsorted!
            false,
        );

        // Now, sort and simplify the tables we got from the sim:
        tables.sort_tables(forrustts::TableSortingFlags::empty());
        let mut samples = SamplesInfo::new();
        for (i, n) in tables.nodes().iter().enumerate() {
            if n.time == num_generations {
                samples.samples.push(i as IdType);
            }
        }

        let mut output = SimplificationOutput::new();
        simplify_tables_without_state(
            &samples,
            SimplificationFlags::empty(),
            &mut tables,
            &mut output,
        )
        .unwrap();

        is_sample = vec![0; tables.num_nodes()];
        for i in is_sample.iter_mut().take(500) {
            *i = 1;
        }

        let simplified_rust_tables = forrustts::tskit_tools::convert_to_tskit_minimal(
            &tables,
            &is_sample,
            forrustts::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        tsk_tables
            .full_sort(tskit::TableSortOptions::default())
            .unwrap();
        tsk_tables
            .simplify(
                &samples.samples,
                tskit::SimplificationOptions::default(),
                false,
            )
            .unwrap();

        // Get tree sequences now
        let tsk_ts = tsk_tables
            .tree_sequence(tskit::TreeSequenceFlags::BUILD_INDEXES)
            .unwrap();

        let rust_ts = simplified_rust_tables
            .tree_sequence(tskit::TreeSequenceFlags::BUILD_INDEXES)
            .unwrap();

        assert_eq!(500, tsk_ts.num_samples());
        assert_eq!(500, rust_ts.num_samples());
        let ne = tsk_ts.edges().num_rows();
        let ne2 = rust_ts.edges().num_rows();
        assert_eq!(ne, ne2);
        let nn = tsk_ts.nodes().num_rows();
        let nn2 = rust_ts.nodes().num_rows();
        assert_eq!(nn, nn2);

        let kc = tsk_ts.kc_distance(&&rust_ts, 0.).unwrap();
        assert!((kc - 0.).abs() < f64::EPSILON);
    }

    #[test]
    fn test_buffer_vs_sort() {
        let num_generations = 5000;
        let genome_length = 1000000;

        let flags = SimulationFlags::USE_STATE;
        let (tables_sorted, is_sample_sorted) = simulate_data(
            num_generations,
            genome_length,
            0.0,
            14613641,
            Some(100),
            flags,
        )
        .unwrap();

        let flags = SimulationFlags::BUFFER_EDGES;
        let (tables_buffered, is_sample_buffered) = simulate_data(
            num_generations,
            genome_length,
            0.0,
            14613641,
            Some(100),
            flags,
        )
        .unwrap();

        // The sums of node times should be the same if we've got
        // both methods working
        let sum_times_sorted: Time = tables_sorted.nodes().iter().map(|x| x.time).sum();
        let sum_times_buffered: Time = tables_buffered.nodes().iter().map(|x| x.time).sum();
        assert_eq!(sum_times_sorted, sum_times_buffered);

        let tables_sorted_tskit = forrustts::tskit_tools::convert_to_tskit_minimal(
            &tables_sorted,
            &is_sample_sorted,
            forrustts::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let tables_buffered_tskit = forrustts::tskit_tools::convert_to_tskit_minimal(
            &tables_buffered,
            &is_sample_buffered,
            forrustts::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let sorted_ts = tables_sorted_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        let buffered_ts = tables_buffered_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        assert_eq!(500, sorted_ts.num_samples());
        assert_eq!(500, buffered_ts.num_samples());
        let kc = sorted_ts.kc_distance(&&buffered_ts, 0.).unwrap();
        assert!((kc - 0.).abs() < f64::EPSILON);
    }

    // The KC distance code will barf on trees where samples
    // have different ages.  So we have to use less direct methods
    // to compare.
    #[test]
    fn test_buffer_vs_sort_overlapping_generations() {
        let num_generations = 5000;
        let genome_length = 1000000;

        let flags = SimulationFlags::USE_STATE;
        let (tables_sorted, is_sample_sorted) = simulate_data(
            num_generations,
            genome_length,
            0.5,
            14613641,
            Some(100),
            flags,
        )
        .unwrap();

        let flags = SimulationFlags::BUFFER_EDGES;
        let (tables_buffered, is_sample_buffered) = simulate_data(
            num_generations,
            genome_length,
            0.5,
            14613641,
            Some(100),
            flags,
        )
        .unwrap();

        assert_eq!(tables_sorted.num_nodes(), tables_buffered.num_nodes());

        // The sums of node times should be the same if we've got
        // both methods working
        let sum_times_sorted: Time = tables_sorted.nodes().iter().map(|x| x.time).sum();
        let sum_times_buffered: Time = tables_buffered.nodes().iter().map(|x| x.time).sum();
        assert_eq!(sum_times_sorted, sum_times_buffered);

        let tables_sorted_tskit = forrustts::tskit_tools::convert_to_tskit_minimal(
            &tables_sorted,
            &is_sample_sorted,
            forrustts::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let tables_buffered_tskit = forrustts::tskit_tools::convert_to_tskit_minimal(
            &tables_buffered,
            &is_sample_buffered,
            forrustts::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let sorted_ts = tables_sorted_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        let buffered_ts = tables_buffered_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        assert_eq!(500, sorted_ts.num_samples());
        assert_eq!(500, buffered_ts.num_samples());
        assert_eq!(sorted_ts.num_trees(), buffered_ts.num_trees());
    }
}
