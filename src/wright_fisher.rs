//! Examples of tree sequence recording with neutral
//! Wright-Fisher models.
//!
//! This module provides a means of generating testing
//! code and benchmarking utilities.  However, some of
//! the concepts here that are *not* public may be useful
//! to others.  Feel free to copy them!
use crate::simplify_from_edge_buffer::simplify_from_edge_buffer;
use crate::simplify_tables::*;
use crate::tables::{validate_edge_table, TableCollection};
use crate::tsdef::*;
use crate::EdgeBuffer;
use crate::ForrusttsError;
use crate::Segment;
use crate::SimplificationBuffers;
use crate::SimplificationFlags;
use crate::SimplificationOutput;
use bitflags::bitflags;
use rgsl::rng::algorithms::mt19937;

// Some of the material below seems like a candidate for a public API,
// but we need to decide here if this package should provide that.
// If so, then many of these types should not be here, as they have nothing
// to do with Wright-Fisher itself, and are instead more general.

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
    p0node0: IdType,
    p0node1: IdType,
    p1node0: IdType,
    p1node1: IdType,
}

impl Birth {
    pub const fn new(index: usize, parent0: &Parent, parent1: &Parent) -> Birth {
        Birth {
            index,
            p0node0: parent0.node0,
            p0node1: parent0.node1,
            p1node0: parent1.node0,
            p1node1: parent1.node1,
        }
    }
}

type VecParent = Vec<Parent>;
type VecBirth = Vec<Birth>;

struct PopulationState {
    pub parents: VecParent,
    pub births: VecBirth,
    pub alive_at_last_simplification: Vec<IdType>,
    pub edge_buffer: EdgeBuffer,
    pub tables: TableCollection,
}

impl PopulationState {
    pub fn new(genome_length: Position) -> Self {
        PopulationState {
            parents: vec![],
            births: vec![],
            alive_at_last_simplification: vec![],
            edge_buffer: EdgeBuffer::new(),
            tables: TableCollection::new(genome_length).unwrap(),
        }
    }
}

fn deaths_and_parents(psurvival: f64, rng: &mut rgsl::Rng, pop: &mut PopulationState) {
    pop.births.clear();
    for i in 0..pop.parents.len() {
        if rng.uniform() > psurvival {
            let parent0 = rng.flat(0., pop.parents.len() as f64) as usize;
            let parent1 = rng.flat(0., pop.parents.len() as f64) as usize;
            pop.births
                .push(Birth::new(i, &pop.parents[parent0], &pop.parents[parent1]));
        }
    }
}

/// Decide which node to pass on from a parent.
fn mendel(rng: &mut rgsl::Rng, n0: IdType, n1: IdType) -> (IdType, IdType) {
    if rng.uniform() < 0.5 {
        (n1, n0)
    } else {
        (n0, n1)
    }
}

fn generate_births(
    littler: f64,
    birth_time: Time,
    rng: &mut rgsl::Rng,
    breakpoints: &mut Vec<Position>,
    pop: &mut PopulationState,
    recorder: impl Fn((IdType, IdType), IdType, &[Position], &mut TableCollection, &mut EdgeBuffer),
) {
    for b in &pop.births {
        let parent0_nodes = mendel(rng, b.p0node0, b.p0node1);
        let parent1_nodes = mendel(rng, b.p1node0, b.p1node1);

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

fn next_breakpoint_distance(
    v: Position,
    breakpoints: &[Position],
    f: impl Fn(Position, Position) -> bool,
) -> usize {
    let i = match breakpoints.iter().position(|x| f(*x, v)) {
        Some(x) => x,
        None => breakpoints.len(),
    };
    i
}

/// If some breakpoints are not unique,
/// prune the input to only include those
/// occurring an odd number of times.
/// This is needed because:
/// Only breakpoints occurring odd numbers of times
/// affect the offspring gamete.
/// We need to ensure we don't do things like
/// add edges with left == right, etc..
fn prune_breakpoints(breakpoints: &mut Vec<Position>) {
    let mut i: usize = 1;
    while i < breakpoints.len() {
        if breakpoints[i - 1] == breakpoints[i] {
            i -= 1;
            break;
        }
        i += 1;
    }

    if i < breakpoints.len() {
        let mut odd_breakpoints = Vec::<Position>::new();
        let mut start: usize = 0;
        while i < breakpoints.len() {
            let not_equal = next_breakpoint_distance(
                breakpoints[i],
                &breakpoints[i..breakpoints.len()],
                |a, b| a != b,
            );
            let even = not_equal % 2 == 0;
            for j in breakpoints.iter().take(i + 1 - even as usize).skip(start) {
                odd_breakpoints.push(*j);
            }
            start = i + not_equal;
            if start >= breakpoints.len() {
                break;
            }
            i = start
                + next_breakpoint_distance(
                    breakpoints[i],
                    &breakpoints[i..breakpoints.len()],
                    |a, b| a == b,
                );
        }
        std::mem::swap(breakpoints, &mut odd_breakpoints);
    }
}

fn recombination_breakpoints(
    littler: f64,
    maxlen: Position,
    rng: &mut rgsl::Rng,
    breakpoints: &mut Vec<Position>,
) {
    breakpoints.clear();
    let nxovers = rng.poisson(littler);
    for _ in 0..nxovers {
        breakpoints.push(rng.flat(0., maxlen as f64) as Position);
    }
    breakpoints.sort_unstable();
    prune_breakpoints(breakpoints);
    if !breakpoints.is_empty() {
        breakpoints.push(Position::MAX);
    }
}

fn fill_samples(parents: &[Parent], samples: &mut Vec<IdType>) {
    samples.clear();
    for p in parents {
        samples.push(p.node0);
        samples.push(p.node1);
    }
}

fn sort_and_simplify(
    flags: SimulationFlags,
    samples: &[IdType],
    state: &mut SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut SimplificationOutput,
) {
    if !flags.contains(SimulationFlags::BUFFER_EDGES) {
        pop.tables.sort_tables_for_simplification();
        debug_assert!(validate_edge_table(
            pop.tables.genome_length(),
            pop.tables.edges(),
            pop.tables.nodes()
        )
        .unwrap());
        if flags.contains(SimulationFlags::USE_STATE) {
            simplify_tables(
                samples,
                SimplificationFlags::empty(),
                state,
                &mut pop.tables,
                output,
            )
            .unwrap();
        } else {
            simplify_tables_without_state(
                samples,
                SimplificationFlags::empty(),
                &mut pop.tables,
                output,
            )
            .unwrap();
        }
        debug_assert!(validate_edge_table(
            pop.tables.genome_length(),
            pop.tables.edges(),
            pop.tables.nodes()
        )
        .unwrap());
    } else {
        simplify_from_edge_buffer(
            samples,
            &pop.alive_at_last_simplification,
            SimplificationFlags::empty(),
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
    samples: &mut Vec<IdType>,
    state: &mut SimplificationBuffers,
    pop: &mut PopulationState,
    output: &mut SimplificationOutput,
) {
    fill_samples(&pop.parents, samples);
    sort_and_simplify(flags, samples, state, pop, output);

    for p in &mut pop.parents {
        p.node0 = output.idmap[p.node0 as usize];
        p.node1 = output.idmap[p.node1 as usize];
    }

    if flags.contains(SimulationFlags::BUFFER_EDGES) {
        pop.alive_at_last_simplification.clear();
        for p in &pop.parents {
            pop.alive_at_last_simplification.push(p.node0);
            pop.alive_at_last_simplification.push(p.node1);
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
    pub seed: usize,
    /// How many birth steps to simulate.
    /// If [``PopulationParams::psurvival``] is 0.0,
    /// then this is the number of generations in the
    /// standard Wright-Fisher model.
    pub nsteps: Time,
    /// Bitwise flag tweaking the behavior of the
    /// simplification algorithm.
    pub flags: SimulationFlags,
}

impl SimulationParams {
    /// Create a new instance
    pub fn new(
        simplification_interval: Option<Time>,
        seed: usize,
        nsteps: Time,
        flags: SimulationFlags,
    ) -> Self {
        SimulationParams {
            simplification_interval,
            seed,
            nsteps,
            flags,
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

    let mut rng;

    match rgsl::Rng::new(mt19937()) {
        None => panic!("failed to allocate rng!"),
        Some(x) => rng = x,
    }

    rng.set(params.seed);

    let mut pop = PopulationState::new(pop_params.genome_length);
    let mut samples: Vec<IdType> = vec![];
    let mut breakpoints = vec![];

    // Record nodes for the first generation
    // Nodes will have birth time 0 in deme 0.
    for i in 0..pop_params.size {
        let n0 = pop.tables.add_node(0, 0).unwrap();
        let n1 = pop.tables.add_node(0, 0).unwrap();
        pop.parents.push(Parent::new(i as usize, n0, n1));
    }

    for i in 0..pop.tables.num_nodes() {
        pop.alive_at_last_simplification.push(i as IdType);
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

    #[test]
    fn test_prune_breakpoints() {
        let mut b = vec![1, 2, 3, 3, 4];
        prune_breakpoints(&mut b);
        assert_eq!(b.len(), 3);
        assert!(b == vec![1, 2, 4]);

        b = vec![1, 1, 2, 3, 3, 4, 4, 5];
        prune_breakpoints(&mut b);
        assert_eq!(b.len(), 2);
        assert!(b == vec![2, 5]);

        b = vec![1, 1, 2, 3, 3, 3, 4, 4, 5];
        prune_breakpoints(&mut b);
        assert_eq!(b.len(), 3);
        assert!(b == vec![2, 3, 5]);
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

        for (i, j) in tables.nodes_.iter().zip(tables_state.nodes_) {
            assert_eq!(i.time, j.time);
            assert_eq!(i.deme, j.deme);
        }
        for (i, j) in tables.edges_.iter().zip(tables_state.edges_) {
            assert_eq!(i.left, j.left);
            assert_eq!(i.right, j.right);
            assert_eq!(i.parent, j.parent);
            assert_eq!(i.child, j.child);
        }
    }
}
