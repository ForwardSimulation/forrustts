use crate::simplify_tables::{simplify_tables, simplify_tables_with_buffers};
use crate::tables::{TableCollection, TreeSequenceRecordingInterface};
use crate::SimplificationBuffers;
use rgsl;
use rgsl::rng::algorithms::mt19937;

// Some of the material below seems like a candidate for a public API,
// but we need to decide here if this package should provide that.
// If so, then many of these types should not be here, as they have nothing
// to do with Wright-Fisher itself, and are instead more general.

//FIXME:
type TsInt = i32;
type SamplesVec = Vec<TsInt>;

struct Parent {
    index: usize,
    node0: TsInt,
    node1: TsInt,
}

impl Parent {
    pub const fn new(index: usize, node0: TsInt, node1: TsInt) -> Parent {
        return Parent {
            index: index,
            node0: node0,
            node1: node1,
        };
    }
}

struct Birth {
    index: usize,
    p0node0: TsInt,
    p0node1: TsInt,
    p1node0: TsInt,
    p1node1: TsInt,
}

impl Birth {
    pub const fn new(index: usize, parent0: &Parent, parent1: &Parent) -> Birth {
        return Birth {
            index: index,
            p0node0: parent0.node0,
            p0node1: parent0.node1,
            p1node0: parent1.node0,
            p1node1: parent1.node1,
        };
    }
}

type VecParent = Vec<Parent>;
type VecBirth = Vec<Birth>;

fn deaths_and_parents(
    parents: &VecParent,
    psurvival: f64,
    rng: &mut rgsl::Rng,
    births: &mut VecBirth,
) -> () {
    births.clear();
    for i in 0..parents.len() {
        if rng.uniform() > psurvival {
            let parent0 = rng.flat(0., parents.len() as f64) as usize;
            let parent1 = rng.flat(0., parents.len() as f64) as usize;
            births.push(Birth::new(i, &parents[parent0], &parents[parent1]));
        }
    }
}

/// Decide which node to pass on from a parent.
fn mendel(rng: &mut rgsl::Rng, n0: TsInt, n1: TsInt) -> (TsInt, TsInt) {
    if rng.uniform() < 0.5 {
        return (n1, n0);
    }
    return (n0, n1);
}

fn generate_births(
    births: &VecBirth,
    littler: f64,
    birth_time: i64,
    rng: &mut rgsl::Rng,
    parents: &mut VecParent,
    breakpoints: &mut Vec<i64>,
    tables: &mut TableCollection,
) -> () {
    for b in births {
        let parent0_nodes = mendel(rng, b.p0node0, b.p0node1);
        let parent1_nodes = mendel(rng, b.p1node0, b.p1node1);

        // Record 2 new nodes
        let new_node_0: TsInt = tables.add_node(birth_time, 0).unwrap();
        let new_node_1: TsInt = tables.add_node(birth_time, 0).unwrap();

        recombination_breakpoints(littler, tables.get_length(), rng, breakpoints);
        record_edges(parent0_nodes, new_node_0, breakpoints, tables);

        recombination_breakpoints(littler, tables.get_length(), rng, breakpoints);
        record_edges(parent1_nodes, new_node_1, breakpoints, tables);

        parents[b.index].index = b.index;
        parents[b.index].node0 = new_node_0;
        parents[b.index].node1 = new_node_1;
    }
}

fn record_edges(
    parents: (TsInt, TsInt),
    child: TsInt,
    breakpoints: &Vec<i64>,
    tables: &mut TableCollection,
) -> () {
    if breakpoints.len() == 0 {
        tables
            .add_edge(0, tables.get_length(), parents.0, child)
            .unwrap();
        return;
    }

    // If we don't have a breakpoint at 0, add an edge
    if breakpoints[0] != 0 {
        tables
            .add_edge(0, breakpoints[0], parents.0, child)
            .unwrap();
    }

    // FIXME: this will generate invalid edge tables when breakpoints
    // are not unique.  Gotta mimic what fwdpp does here.
    for i in 1..breakpoints.len() {
        let a = breakpoints[i - 1];
        let b = if i < (breakpoints.len() - 1) {
            breakpoints[i]
        } else {
            tables.get_length()
        };
        if i % 2 == 0 {
            tables.add_edge(a, b, parents.0, child).unwrap();
        } else {
            tables.add_edge(a, b, parents.1, child).unwrap();
        }
    }
}

fn next_breakpoint_distance(v: i64, breakpoints: &[i64], f: impl Fn(i64, i64) -> bool) -> usize {
    let i = match breakpoints.iter().position(|x| f(*x, v)) {
        Some(x) => x,
        None => breakpoints.len(),
    };
    return i;
}

/// If some breakpoints are not unique,
/// prune the input to only include those
/// occurring an odd number of times.
/// This is needed because:
/// Only breakpoints occurring odd numbers of times
/// affect the offspring gamete.
/// We need to ensure we don't do things like
/// add edges with left == right, etc..
fn prune_breakpoints(breakpoints: &mut Vec<i64>) {
    let mut i: usize = 1;
    while i < breakpoints.len() {
        if breakpoints[i - 1] == breakpoints[i] {
            i = i - 1;
            break;
        }
        i += 1;
    }

    if i < breakpoints.len() {
        let mut odd_breakpoints = Vec::<i64>::new();
        let mut start: usize = 0;
        while i < breakpoints.len() {
            let not_equal = next_breakpoint_distance(
                breakpoints[i],
                &breakpoints[i..breakpoints.len()],
                |a, b| {
                    return a != b;
                },
            );
            let even = if not_equal % 2 == 0 { true } else { false };
            for j in start..i + 1 - even as usize {
                odd_breakpoints.push(breakpoints[j]);
            }
            start = i + not_equal;
            if start >= breakpoints.len() {
                break;
            }
            i = start
                + next_breakpoint_distance(
                    breakpoints[i],
                    &breakpoints[i..breakpoints.len()],
                    |a, b| {
                        return a == b;
                    },
                );
        }
        std::mem::swap(breakpoints, &mut odd_breakpoints);
    }
}

fn recombination_breakpoints(
    littler: f64,
    maxlen: i64,
    rng: &mut rgsl::Rng,
    breakpoints: &mut Vec<i64>,
) -> () {
    breakpoints.clear();
    let nxovers = rng.poisson(littler);
    for _ in 0..nxovers {
        breakpoints.push(rng.flat(0., maxlen as f64) as i64);
    }
    breakpoints.sort();
    prune_breakpoints(breakpoints);
    if breakpoints.len() > 0 {
        breakpoints.push(std::i64::MAX);
    }
}

// NOTE: I've apparently decided on changing my naming convention?
fn fill_samples(parents: &VecParent, samples: &mut SamplesVec) -> () {
    samples.clear();
    for p in parents {
        samples.push(p.node0);
        samples.push(p.node1);
    }
}

fn sort_and_simplify(
    use_state: bool,
    samples: &SamplesVec,
    state: &mut SimplificationBuffers,
    tables: &mut TableCollection,
) -> SamplesVec {
    tables.sort_tables_for_simplification();
    debug_assert!(
        tables.validate_edge_table().unwrap()
    );
    let idmap = if use_state == true {
        simplify_tables_with_buffers(samples, state, tables)
    } else {
        simplify_tables(samples, tables)
    };
    debug_assert!(
        tables.validate_edge_table().unwrap()
    );
    return idmap;
}

fn simplify_and_remap_nodes(
    use_state: bool,
    samples: &mut SamplesVec,
    parents: &mut VecParent,
    state: &mut SimplificationBuffers,
    tables: &mut TableCollection,
) -> () {
    fill_samples(parents, samples);
    let idmap = sort_and_simplify(use_state, samples, state, tables);
    for p in parents {
        p.node0 = idmap[p.node0 as usize];
        p.node1 = idmap[p.node1 as usize];
    }
}

fn validate_simplification_interval(x: i64) -> i64 {
    if x < 1 {
        panic!("simplification_interval must be None or >= 1");
    }
    return x;
}

// NOTE: this function is a copy of the simulation
// found in fwdpp/examples/edge_buffering.cc

fn neutral_wf_impl(
    seed: usize,
    popsize: u32,
    nsteps: i64,
    genome_length: i64,
    littler: f64,
    psurvival: f64,
    simplification_interval: Option<i64>,
    use_state: bool,
) -> TableCollection {
    // FIXME: gotta validate input params!

    let mut actual_simplification_interval: i64 = -1;

    match simplification_interval {
        None => (),
        Some(x) => actual_simplification_interval = validate_simplification_interval(x),
    }

    let mut rng;

    match rgsl::Rng::new(mt19937()) {
        None => panic!("failed to allocate rng!"),
        Some(x) => rng = x,
    }

    rng.set(seed);

    let mut tables = TableCollection::new(genome_length).unwrap();
    let mut parents = VecParent::new();
    let mut births = VecBirth::new();
    let mut samples = SamplesVec::new();
    let mut breakpoints = vec![];

    // Record nodes for the first generation
    // Nodes will have birth time 0 in deme 0.
    for i in 0..popsize {
        let n0 = tables.add_node(0, 0).unwrap();
        let n1 = tables.add_node(0, 0).unwrap();
        parents.push(Parent::new(i as usize, n0, n1));
    }

    let mut simplified = false;
    let mut state = SimplificationBuffers::new();

    for birth_time in 1..(nsteps + 1) {
        deaths_and_parents(&parents, psurvival, &mut rng, &mut births);
        generate_births(
            &births,
            littler,
            birth_time,
            &mut rng,
            &mut parents,
            &mut breakpoints,
            &mut tables,
        );
        if actual_simplification_interval != -1 && birth_time % actual_simplification_interval == 0
        {
            simplify_and_remap_nodes(
                use_state,
                &mut samples,
                &mut parents,
                &mut state,
                &mut tables,
            );
            simplified = true;
        } else {
            simplified = false;
        }
    }

    if simplified == false && actual_simplification_interval != -1 {
        simplify_and_remap_nodes(
            use_state,
            &mut samples,
            &mut parents,
            &mut state,
            &mut tables,
        );
    }

    return tables;
}

pub fn neutral_wf(
    seed: usize,
    popsize: u32,
    nsteps: i64,
    genome_length: i64,
    littler: f64,
    psurvival: f64,
    simplification_interval: Option<i64>,
) -> TableCollection {
    return neutral_wf_impl(
        seed,
        popsize,
        nsteps,
        genome_length,
        littler,
        psurvival,
        simplification_interval,
        true,
    );
}

#[cfg(test)]
mod test {

    use super::*;

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

    fn run_sim(use_state: bool) -> TableCollection {
        return neutral_wf_impl(666, 1000, 2000, 100000, 5e-3, 0.0, Some(100), use_state);
    }

    #[test]
    fn compare_state_to_no_state() {
        let tables = run_sim(false);
        let tables_state = run_sim(true);

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
