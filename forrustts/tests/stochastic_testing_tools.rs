#[path = "./neutral_wright_fisher.rs"]
mod neutral_wright_fisher;

use forrustts::*;
pub use neutral_wright_fisher::*;
use std::convert::TryInto;

pub struct SimulatorIterator {
    rng: Rng,
    params: SimulationParams,
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
            rng: Rng::new(params.seed),
            params,
            nreps,
            rep: 0,
        }
    }
}

impl Iterator for SimulatorIterator {
    type Item = SimResults;

    fn next(&mut self) -> Option<Self::Item> {
        if self.rep < self.nreps {
            let seed = self.rng.0.flat(0., usize::MAX as f64) as usize;
            let mut params = self.params;
            params.seed = seed;
            let (mut tables, is_sample) = neutral_wf(params).unwrap();

            tables.sort_tables(TableSortingFlags::empty());

            let mut tsk_tables = forrustts_tskit::export_tables(
                tables.clone(),
                forrustts_tskit::simple_time_reverser(params.nsteps),
                // Do not index tables here!
                // Things are unsorted!
                None,
            )
            .unwrap();
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
