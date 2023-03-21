use forrustts_core::newtypes::Position;
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
    poisson_regions: PoissonRegions,
}

impl GeneticMap {
    pub fn new_from_builder(builder: GeneticMapBuilder) -> Option<Self> {
        let poisson_regions = PoissonRegions::new(builder.poisson())?;
        Some(Self { poisson_regions })
    }
}

impl<T> GenerateBreakpoints<T> for GeneticMap
where
    T: Rng,
{
    fn generate_breakpoints(&mut self, rng: &mut T) {
        unimplemented!("not yet!");
    }
    fn breakpoints(&self) -> &[forrustts_genetic_maps::Breakpoint] {
        unimplemented!("not yet");
    }
}
