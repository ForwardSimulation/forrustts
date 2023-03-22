use std::collections::binary_heap;

use forrustts_core::newtypes::Position;
use forrustts_genetic_maps::BernoulliCrossover;
use forrustts_genetic_maps::Breakpoint;
use forrustts_genetic_maps::GenerateBreakpoints;
use forrustts_genetic_maps::GeneticMapBuilder;
use forrustts_genetic_maps::IndependentAssortment;
use forrustts_genetic_maps::PoissonCrossover;

use rand::Rng;

#[derive(Debug)]
struct PoissonRegions {
    regions: Vec<rand_distr::Uniform<i64>>,
    lookup: rand_distr::WeightedAliasIndex<f64>,
    poisson: rand_distr::Poisson<f64>,
}

impl PoissonRegions {
    fn new(poisson: &[PoissonCrossover]) -> Option<Self> {
        let mut regions = vec![];
        let mut weights = vec![];
        let mut sum_poisson_means = 0.0;
        for p in poisson {
            let mean = p.mean();
            if mean > 0.0 {
                let u = rand_distr::Uniform::new(i64::from(p.left()), i64::from(p.right()));
                regions.push(u);
                weights.push(mean);
                sum_poisson_means += mean;
            }
        }
        let lookup = rand_distr::WeightedAliasIndex::new(weights).ok()?;
        let poisson = rand_distr::Poisson::new(sum_poisson_means).ok()?;
        Some(Self {
            regions,
            lookup,
            poisson,
        })
    }
}

#[derive(Debug)]
struct BernoulliRegions {
    regions: Vec<rand_distr::Uniform<i64>>,
    probabilities: Vec<rand_distr::Bernoulli>,
}

impl BernoulliRegions {
    fn new(binomial: &[BernoulliCrossover]) -> Option<Self> {
        let mut regions = vec![];
        let mut probabilities = vec![];

        for b in binomial {
            let prob = b.probability();
            if prob > 0.0 {
                let u = rand_distr::Uniform::new(i64::from(b.left()), i64::from(b.right()));
                regions.push(u);
                let dist = rand_distr::Bernoulli::new(prob).ok()?;
                probabilities.push(dist);
            }
        }

        Some(Self {
            regions,
            probabilities,
        })
    }
}

#[derive(Debug)]
pub struct GeneticMap {
    poisson_regions: Option<PoissonRegions>,
    bernoulli_regions: Option<BernoulliRegions>,
    independent_assortment: Vec<IndependentAssortment>,
    breakpoints: Vec<Breakpoint>,
}

impl GeneticMap {
    pub fn new_from_builder(builder: GeneticMapBuilder) -> Option<Self> {
        let poisson_regions = PoissonRegions::new(builder.poisson());
        let bernoulli_regions = BernoulliRegions::new(builder.binomial());
        let independent_assortment = builder.independent_assortment;
        Some(Self {
            poisson_regions,
            bernoulli_regions,
            independent_assortment,
            breakpoints: vec![],
        })
    }
}

impl<T> GenerateBreakpoints<T> for GeneticMap
where
    T: Rng,
{
    fn generate_breakpoints(&mut self, rng: &mut T) {
        self.breakpoints.clear();
        if let Some(poisson) = self.poisson_regions.as_ref() {
            let num: u32 = rng.sample(poisson.poisson) as u32;
            for _ in 0..num {
                let idx = rng.sample(&poisson.lookup);
                let pos = rng.sample(poisson.regions[idx]);
                self.breakpoints.push(Breakpoint::Crossover(pos.into()));
            }
        }
        if let Some(bernoulli) = self.bernoulli_regions.as_ref() {
            for (b, r) in bernoulli.probabilities.iter().zip(bernoulli.regions.iter()) {
                let x = rng.sample(b);
                if x {
                    let pos = rng.sample(r);
                    self.breakpoints.push(Breakpoint::Crossover(pos.into()));
                }
            }
        }
        let uint = rand_distr::Uniform::new(0., 1.);
        for i in self.independent_assortment.iter() {
            if rng.sample(uint) <= 0.5 {
                self.breakpoints.push(i.at());
            }
        }
        self.breakpoints.sort();
    }

    fn breakpoints(&self) -> &[forrustts_genetic_maps::Breakpoint] {
        &self.breakpoints
    }
}
