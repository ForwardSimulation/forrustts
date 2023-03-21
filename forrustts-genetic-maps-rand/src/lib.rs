use forrustts_core::newtypes::Position;
use forrustts_genetic_maps::Breakpoint;
use forrustts_genetic_maps::GenerateBreakpoints;
use forrustts_genetic_maps::GeneticMapBuilder;
use forrustts_genetic_maps::PoissonCrossover;

use rand::Rng;

#[derive(Default, Debug)]
struct Regions {
    left: Vec<Position>,
    right: Vec<Position>,
}

#[derive(Debug)]
struct PoissonRegions {
    regions: Regions,
    lookup: rand_distr::WeightedAliasIndex<f64>,
    sum_poisson_means: f64,
}

impl PoissonRegions {
    fn new(poisson: &[PoissonCrossover]) -> Option<Self> {
        let mut regions = Regions::default();
        let mut weights = vec![];
        let mut sum_poisson_means = 0.0;
        for p in poisson {
            let mean = p.mean();
            if mean > 0.0 {
                regions.left.push(p.left());
                regions.right.push(p.right());
                weights.push(mean);
                sum_poisson_means += mean;
            }
        }
        let lookup = rand_distr::WeightedAliasIndex::new(weights).ok()?;
        Some(Self {
            regions,
            lookup,
            sum_poisson_means,
        })
    }
}

#[derive(Debug)]
pub struct GeneticMap {
    poisson_regions: Option<PoissonRegions>,
    breakpoints: Vec<Breakpoint>,
}

impl GeneticMap {
    pub fn new_from_builder(builder: GeneticMapBuilder) -> Option<Self> {
        let poisson_regions = PoissonRegions::new(builder.poisson());
        Some(Self {
            poisson_regions,
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
            // NOTE: all of this object construction should be in PoissonRegions...
            let p = rand_distr::Poisson::new(poisson.sum_poisson_means).unwrap();
            let num: u32 = rng.sample(p) as u32;
            for _ in 0..num {
                let idx = rng.sample(&poisson.lookup);
                let u = rand_distr::Uniform::new(
                    i64::from(poisson.regions.left[idx]),
                    i64::from(poisson.regions.right[idx]),
                );
                let pos = rng.sample(u);
                self.breakpoints.push(Breakpoint::Crossover(pos.into()));
            }
        }
        self.breakpoints.sort();
    }

    fn breakpoints(&self) -> &[forrustts_genetic_maps::Breakpoint] {
        &self.breakpoints
    }
}
