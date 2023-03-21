use forrustts_core::newtypes::Position;

use rand::Rng;

pub trait CrossoverRegionModel {}

pub struct BinomialRegionModel {}
pub struct PoissonRegionModel {}

impl CrossoverRegionModel for BinomialRegionModel {}
impl CrossoverRegionModel for PoissonRegionModel {}

pub struct CrossoverRegion<T>
where
    T: CrossoverRegionModel,
{
    left: Position,
    right: Position,
    parameter: f64,
    marker: std::marker::PhantomData<T>,
}

impl CrossoverRegion<PoissonRegionModel> {
    pub fn new(left: Position, right: Position, mean: f64) -> Option<Self> {
        if left < 0 || right < 0 || right <= left || mean < 0.0 {
            None
        } else {
            Some(Self {
                left,
                right,
                parameter: mean,
                marker: std::marker::PhantomData,
            })
        }
    }
}

pub struct GeneticMap {}

