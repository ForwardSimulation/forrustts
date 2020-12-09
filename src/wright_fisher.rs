use crate::simplify_tables::simplify_tables;
use crate::tables::{validate_edge_table, TableCollection, TablesError};
use crate::tsdef::*;
use rgsl;
use rgsl::rng::algorithms::mt19937;

// Some of the material below seems like a candidate for a public API,
// but we need to decide here if this package should provide that.
// If so, then many of these types should not be here, as they have nothing
// to do with Wright-Fisher itself, and are instead more general.

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

// FIXME: need to deal w/recombination!
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

fn sort_and_simplify(samples: &SamplesVec, tables: &mut TableCollection) -> SamplesVec {
    tables.sort_tables_for_simplification();
    debug_assert!(
        validate_edge_table(tables.get_length(), tables.edges(), tables.nodes()).unwrap()
    );
    let idmap = simplify_tables(samples, tables);
    debug_assert!(
        validate_edge_table(tables.get_length(), tables.edges(), tables.nodes()).unwrap()
    );
    return idmap;
}

fn simplify_and_remap_nodes(
    samples: &mut SamplesVec,
    parents: &mut VecParent,
    tables: &mut TableCollection,
) -> () {
    fill_samples(parents, samples);
    let idmap = sort_and_simplify(samples, tables);
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

pub fn neutral_wf(
    seed: usize,
    popsize: u32,
    nsteps: i64,
    genome_length: i64,
    littler: f64,
    psurvival: f64,
    simplification_interval: Option<i64>,
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
            simplify_and_remap_nodes(&mut samples, &mut parents, &mut tables);
            simplified = true;
        } else {
            simplified = false;
        }
    }

    if simplified == false && actual_simplification_interval != -1 {
        simplify_and_remap_nodes(&mut samples, &mut parents, &mut tables);
    }

    return tables;
}
