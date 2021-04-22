//! Examples of tree sequence recording with neutral
//! Wright-Fisher models.
//!
//! This module provides a means of generating testing
//! code and benchmarking utilities.  However, some of
//! the concepts here that are *not* public may be useful
//! to others.  Feel free to copy them!
use bitflags::bitflags;
use clap::{value_t, value_t_or_exit, App, Arg};
use forrustts::simplify_from_edge_buffer;
use forrustts::simplify_tables;
use forrustts::simplify_tables_without_state;
use forrustts::EdgeBuffer;
use forrustts::ForrusttsError;
use forrustts::IdType;
use forrustts::Position;
use forrustts::SamplesInfo;
use forrustts::Segment;
use forrustts::SimplificationBuffers;
use forrustts::SimplificationFlags;
use forrustts::SimplificationOutput;
use forrustts::TableCollection;
use forrustts::Time;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Geometric, Uniform};

// Some of the material below seems like a candidate for a public API,
// but we need to decide here if this package should provide that.
// If so, then many of these types should not be here, as they have nothing
// to do with Wright-Fisher itself, and are instead more general.

#[derive(Copy, Clone)]
struct Parent {
    index: usize,
    node0: IdType,
    node1: IdType,
}

impl Parent {
    pub const fn new(index: usize, node0: IdType, node1: IdType) -> Parent {
        Parent {
            index,
            node0,
            node1,
        }
    }
}

struct Birth {
    index: usize,
    parent0: Parent,
    parent1: Parent,
}

impl Birth {
    pub const fn new(index: usize, parent0: Parent, parent1: Parent) -> Birth {
        Birth {
            index,
            parent0,
            parent1,
        }
    }
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
                pop.births.push(Birth::new(i, parent0, parent1));
            }
            Some(_) => (),
            None => (),
        }
    }
}

/// Decide which node to pass on from a parent.
fn mendel(rng: &mut StdRng, n0: IdType, n1: IdType) -> (IdType, IdType) {
    let x: f64 = rng.gen();
    match x.partial_cmp(&0.5) {
        Some(std::cmp::Ordering::Less) => (n1, n0),
        Some(_) => (n0, n1),
        None => panic!("Unexpected None"),
    }
}

fn generate_births(
    littler: f64,
    birth_time: Time,
    rng: &mut StdRng,
    breakpoints: &mut Vec<Position>,
    pop: &mut PopulationState,
    recorder: impl Fn((IdType, IdType), IdType, &[Position], &mut TableCollection, &mut EdgeBuffer),
) {
    for b in &pop.births {
        let parent0_nodes = mendel(rng, b.parent0.node0, b.parent0.node1);
        let parent1_nodes = mendel(rng, b.parent1.node0, b.parent1.node1);

        // Record 2 new nodes
        let new_node_0: IdType = pop.tables.add_node(birth_time, 0).unwrap();
        let new_node_1: IdType = pop.tables.add_node(birth_time, 0).unwrap();

        recombination_breakpoints(littler, pop.tables.genome_length(), rng, breakpoints);
        recorder(
            parent0_nodes,
            new_node_0,
            breakpoints,
            &mut pop.tables,
            &mut pop.edge_buffer,
        );

        recombination_breakpoints(littler, pop.tables.genome_length(), rng, breakpoints);
        recorder(
            parent1_nodes,
            new_node_1,
            breakpoints,
            &mut pop.tables,
            &mut pop.edge_buffer,
        );

        pop.parents[b.index].index = b.index;
        pop.parents[b.index].node0 = new_node_0;
        pop.parents[b.index].node1 = new_node_1;
    }
}

fn buffer_edges(
    parents: (IdType, IdType),
    child: IdType,
    breakpoints: &[Position],
    tables: &mut TableCollection,
    buffer: &mut EdgeBuffer,
) {
    if breakpoints.is_empty() {
        buffer
            .extend(parents.0, Segment::new(0, tables.genome_length(), child))
            .unwrap();
        return;
    }

    // If we don't have a breakpoint at 0, add an edge
    if breakpoints[0] != 0 {
        buffer
            .extend(parents.0, Segment::new(0, breakpoints[0], child))
            .unwrap();
    }

    for i in 1..breakpoints.len() {
        let a = breakpoints[i - 1];
        let b = if i < (breakpoints.len() - 1) {
            breakpoints[i]
        } else {
            tables.genome_length()
        };
        if i % 2 == 0 {
            buffer.extend(parents.0, Segment::new(a, b, child)).unwrap();
        } else {
            buffer.extend(parents.1, Segment::new(a, b, child)).unwrap();
        }
    }
}

fn record_edges(
    parents: (IdType, IdType),
    child: IdType,
    breakpoints: &[Position],
    tables: &mut TableCollection,
    _: &mut EdgeBuffer,
) {
    if breakpoints.is_empty() {
        tables
            .add_edge(0, tables.genome_length(), parents.0, child)
            .unwrap();
        return;
    }

    // If we don't have a breakpoint at 0, add an edge
    if breakpoints[0] != 0 {
        tables
            .add_edge(0, breakpoints[0], parents.0, child)
            .unwrap();
    }

    for i in 1..breakpoints.len() {
        let a = breakpoints[i - 1];
        let b = if i < (breakpoints.len() - 1) {
            breakpoints[i]
        } else {
            tables.genome_length()
        };
        if i % 2 == 0 {
            tables.add_edge(a, b, parents.0, child).unwrap();
        } else {
            tables.add_edge(a, b, parents.1, child).unwrap();
        }
    }
}

fn recombination_breakpoints(
    littler: f64,
    maxlen: Position,
    rng: &mut StdRng,
    breakpoints: &mut Vec<Position>,
) {
    breakpoints.clear();
    match littler.partial_cmp(&0.0) {
        Some(std::cmp::Ordering::Greater) => {
            let geom = match Geometric::new(littler / maxlen as f64) {
                Ok(g) => g,
                Err(e) => panic!("{}", e),
            };
            let mut current_pos: Position = 0;
            loop {
                let next_length = rng.sample(geom) as Position;
                if current_pos + next_length < maxlen {
                    current_pos += next_length;
                    breakpoints.push(current_pos);
                } else {
                    break;
                }
            }
        }
        Some(_) => {}
        None => (),
    }
    if !breakpoints.is_empty() {
        breakpoints.push(Position::MAX);
    }
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
        pop.tables
            .sort_tables(forrustts::TableSortingFlags::empty());
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

// NOTE: this function is a copy of the simulation
// found in fwdpp/examples/edge_buffering.cc

/// Parameters of a population to be evolved by
/// [``neutral_wf``].
pub struct PopulationParams {
    /// Diploid population size
    pub size: u32,
    /// Genome length.  Easiest to think "base pairs",
    /// but more abstract concepts are valid.
    pub genome_length: Position,
    /// Mean number of crossovers (per mating).
    pub littler: f64,
    /// Survival probability.  Must be 0 <= p < 1.
    pub psurvival: f64,
}

impl PopulationParams {
    /// Create a new instance
    ///
    /// # Limitations
    ///
    /// As of 0.1.0, input values are not validated.
    pub fn new(size: u32, genome_length: Position, littler: f64, psurvival: f64) -> Self {
        PopulationParams {
            size,
            genome_length,
            littler,
            psurvival,
        }
    }
}

bitflags! {
    /// Bitwise flag tweaking the behavior of the
    /// simplification algorithm.
    #[derive(Default)]
    pub struct SimulationFlags: u32
    {
        /// If set, and [``BUFFER_EDGES``] is not set,
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

/// Parameters of a simulation to be executed
/// by [``neutral_wf``].
pub struct SimulationParams {
    /// How often to apply the simplification algorithm.
    /// If ``None``, then simplification never happens.
    /// If the value is ``Some(Time)``, then simplification
    /// will occur after that many time steps.
    pub simplification_interval: Option<Time>,
    /// Random number seed
    pub seed: u64,
    /// How many birth steps to simulate.
    /// If [``PopulationParams::psurvival``] is 0.0,
    /// then this is the number of generations in the
    /// standard Wright-Fisher model.
    pub nsteps: Time,
    /// Bitwise flag tweaking the behavior of the
    /// simplification algorithm.
    pub flags: SimulationFlags,
    /// Flags to affect simplification behavior
    pub simplification_flags: SimplificationFlags,
}

impl SimulationParams {
    /// Create a new instance
    /// [``SimulationParams::simplification_flags``]
    /// is initialized to empty.
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
            simplification_flags: SimplificationFlags::empty(),
        }
    }
}

/// Run a simulation of an idealized population.
///
/// This function simulates a constant-sized population.
/// This is a Wright-Fisher population with the possibility
/// of overlapping generations.
///
/// # Parameters
///
/// * pop_params is an instance of [``PopulationParams``].
/// * params is an instance of [``SimulationParams``].
///
/// # Return values
///
/// The return value is a tuple containing a [``TableCollection``]
/// and a vector.  The length of the vector is equal to the
/// length of the [``crate::NodeTable``] in the return value. The vector
/// contains 1 if the node at that index is alive at the end of
/// the simlation and 0 otherwise.
///
/// # Error
///
/// The simulation may return [``ForrusttsError``] if an error
/// is encountered.
///
/// # Limitations
///
/// As of version 0.1.0, the input values are not
/// thoroughly validated.
pub fn neutral_wf(
    pop_params: PopulationParams,
    params: SimulationParams,
) -> Result<(TableCollection, Vec<i32>), ForrusttsError> {
    // FIXME: gotta validate input params!

    let mut actual_simplification_interval: Time = -1;

    match params.simplification_interval {
        None => (),
        Some(x) => actual_simplification_interval = validate_simplification_interval(x),
    }

    let mut rng = StdRng::seed_from_u64(params.seed);

    let mut pop = PopulationState::new(pop_params.genome_length);
    let mut samples: SamplesInfo = Default::default();
    let mut breakpoints = vec![];

    // Record nodes for the first generation
    // Nodes will have birth time 0 in deme 0.
    for i in 0..pop_params.size {
        let n0 = pop.tables.add_node(0, 0).unwrap();
        let n1 = pop.tables.add_node(0, 0).unwrap();
        pop.parents.push(Parent::new(i as usize, n0, n1));
    }

    for i in 0..pop.tables.num_nodes() {
        samples.edge_buffer_founder_nodes.push(i as IdType);
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
        deaths_and_parents(pop_params.psurvival, &mut rng, &mut pop);
        generate_births(
            pop_params.littler,
            birth_time,
            &mut rng,
            &mut breakpoints,
            &mut pop,
            new_edge_handler,
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
            Arg::with_name("rho")
                .long("rho")
                .help("Scaled rate of crossing over, 4Nr")
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
    let rho = value_t_or_exit!(matches.value_of("rho"), f64);
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

    let r = rho / (4.0 * popsize as f64);

    let flags = if matches.is_present("buffer_edges") {
        SimulationFlags::BUFFER_EDGES
    } else {
        SimulationFlags::USE_STATE
    };

    let mut simplification_flags = SimplificationFlags::empty();

    if validate_tables {
        simplification_flags |= SimplificationFlags::VALIDATE_ALL;
    }

    let (mut tables, is_sample) = neutral_wf(
        PopulationParams {
            size: popsize,
            genome_length: 10000000,
            littler: r,
            psurvival,
        },
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

    fn run_sim(use_state: bool) -> (TableCollection, Vec<i32>) {
        let flags = if use_state {
            SimulationFlags::USE_STATE
        } else {
            SimulationFlags::empty()
        };
        neutral_wf(
            PopulationParams {
                size: 1000,
                genome_length: 100000,
                littler: 5e-3,
                psurvival: 0.0,
            },
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
            PopulationParams {
                size: 250,
                genome_length,
                littler: 5e-3,
                psurvival,
            },
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
