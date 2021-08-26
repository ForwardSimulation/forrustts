#[path = "./neutral_wright_fisher.rs"]
mod neutral_wright_fisher;

use forrustts_tables_trees::*;
pub use neutral_wright_fisher::*;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::Uniform;
use std::convert::TryInto;

pub struct SimulatorIterator {
    rng: StdRng,
    params: SimulationParams,
    make_seed: Uniform<u64>,
    nreps: i32,
    rep: i32,
}

pub struct SimResults {
    pub tables: TableCollection,
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
            let (mut tables, is_sample) = neutral_wf(params).unwrap();

            tables.sort_tables(TableSortingFlags::empty());

            let mut tsk_tables = forrustts_tables_trees::tskit_tools::convert_to_tskit_minimal(
                &tables,
                &is_sample,
                forrustts_tables_trees::tskit_tools::simple_time_reverser(params.nsteps.into()),
                // Do not index tables here!
                // Things are unsorted!
                false,
            );
            add_tskit_mutation_site_tables(&tables, params.nsteps.into(), &mut tsk_tables);
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

pub fn make_samples(l: &[i32]) -> SamplesInfo {
    let mut rv = SamplesInfo::default();
    for (i, j) in l.iter().enumerate() {
        if *j == 1 {
            rv.samples.push(i.try_into().unwrap());
        }
    }
    rv
}
