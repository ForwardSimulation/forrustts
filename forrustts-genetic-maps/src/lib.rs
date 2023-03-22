//! Core types for genetic maps

use forrustts_core::newtypes::Position;
use rand::Rng;

/// Breakpoint positions from crossover events
///
/// # Notes
///
/// * Comparision operations are based
///   on the stored
///   [`Position`](forrustts_core::newtypes::Position).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Breakpoint {
    Crossover(Position),
    IndependentAssortment(Position),
}

impl From<Breakpoint> for Position {
    fn from(value: Breakpoint) -> Self {
        match value {
            Breakpoint::IndependentAssortment(p) => p,
            Breakpoint::Crossover(p) => p,
        }
    }
}

impl PartialOrd for Breakpoint {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let left = Position::from(*self);
        let right = Position::from(*other);
        left.partial_cmp(&right)
    }
}

impl Ord for Breakpoint {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub trait GenerateBreakpoints {
    fn generate_breakpoints<T: Rng>(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Breakpoint];
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PoissonCrossover {
    left: Position,
    right: Position,
    mean: f64,
}

impl PoissonCrossover {
    pub fn new<L, R>(left: L, right: R, mean: f64) -> Option<Self>
    where
        L: Into<Position>,
        R: Into<Position>,
    {
        let left = left.into();
        let right = right.into();
        if left < 0 || right < 0 || right <= left || !mean.is_finite() || mean < 0.0 {
            None
        } else {
            Some(Self { left, right, mean })
        }
    }

    pub fn left(&self) -> Position {
        self.left
    }
    pub fn right(&self) -> Position {
        self.right
    }
    pub fn mean(&self) -> f64 {
        self.mean
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BernoulliCrossover {
    left: Position,
    right: Position,
    probability: f64,
}

impl BernoulliCrossover {
    pub fn new<L, R>(left: L, right: R, probability: f64) -> Option<Self>
    where
        L: Into<Position>,
        R: Into<Position>,
    {
        let left = left.into();
        let right = right.into();
        if left < 0
            || right < 0
            || right <= left
            || !probability.is_finite()
            || !(0.0..=1.0).contains(&probability)
        {
            None
        } else {
            Some(Self {
                left,
                right,
                probability,
            })
        }
    }
    pub fn left(&self) -> Position {
        self.left
    }
    pub fn right(&self) -> Position {
        self.right
    }
    pub fn probability(&self) -> f64 {
        self.probability
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct IndependentAssortment {
    at: Position,
}

impl IndependentAssortment {
    pub fn new<P: Into<Position>>(at: P) -> Option<Self> {
        let at = at.into();
        if at >= 0 {
            Some(Self { at })
        } else {
            None
        }
    }

    pub fn at(&self) -> Breakpoint {
        Breakpoint::IndependentAssortment(self.at)
    }
}

#[derive(Default, Debug, Clone)]
pub struct GeneticMapBuilder {
    poisson: Vec<PoissonCrossover>,
    bernoulli: Vec<BernoulliCrossover>,
    independent_assortment: Vec<IndependentAssortment>,
}

impl GeneticMapBuilder {
    pub fn extend_poisson(mut self, intervals: &[PoissonCrossover]) -> Self {
        self.poisson.extend_from_slice(intervals);
        self
    }

    pub fn extend_bernoulli(mut self, intervals: &[BernoulliCrossover]) -> Self {
        self.bernoulli.extend_from_slice(intervals);
        self
    }

    pub fn extend_independent_assortment(mut self, at: &[IndependentAssortment]) -> Self {
        self.independent_assortment.extend_from_slice(at);
        self
    }

    pub fn poisson(&self) -> &[PoissonCrossover] {
        &self.poisson
    }

    pub fn bernoulli(&self) -> &[BernoulliCrossover] {
        &self.bernoulli
    }

    pub fn independent_assortment(&self) -> &[IndependentAssortment] {
        &self.independent_assortment
    }
}

#[derive(Debug)]
struct PoissonRegions {
    regions: Vec<rand_distr::Uniform<i64>>,
    lookup: rand_distr::WeightedAliasIndex<f64>,
    poisson: rand_distr::Poisson<f64>,
}

impl PoissonRegions {
    fn new(poisson: Vec<PoissonCrossover>) -> Option<Self> {
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

    fn generate<T: Rng>(&self, rng: &mut T) -> Breakpoint {
        let i = rng.sample(&self.lookup);
        Breakpoint::Crossover(rng.sample(self.regions[i]).into())
    }
}

#[derive(Debug)]
struct BernoulliRegions {
    regions: Vec<rand_distr::Uniform<i64>>,
    probabilities: Vec<rand_distr::Bernoulli>,
}

impl BernoulliRegions {
    fn new(binomial: Vec<BernoulliCrossover>) -> Option<Self> {
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

    fn len(&self) -> usize {
        self.regions.len()
    }

    fn generate<T: Rng>(&self, index: usize, rng: &mut T) -> Option<Breakpoint> {
        if rng.sample(self.probabilities[index]) {
            Some(Breakpoint::Crossover(
                rng.sample(self.regions[index]).into(),
            ))
        } else {
            None
        }
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
        let (poisson, bernoulli, independent_assortment) = (
            builder.poisson,
            builder.bernoulli,
            builder.independent_assortment,
        );
        let poisson_regions = PoissonRegions::new(poisson);
        let bernoulli_regions = BernoulliRegions::new(bernoulli);
        Some(Self {
            poisson_regions,
            bernoulli_regions,
            independent_assortment,
            breakpoints: vec![],
        })
    }
}

impl GenerateBreakpoints for GeneticMap {
    fn generate_breakpoints<T: Rng>(&mut self, rng: &mut T) {
        self.breakpoints.clear();
        if let Some(poisson) = self.poisson_regions.as_ref() {
            let num: u32 = rng.sample(poisson.poisson) as u32;
            for _ in 0..num {
                let pos = poisson.generate(rng);
                self.breakpoints.push(pos);
            }
        }
        if let Some(bernoulli) = self.bernoulli_regions.as_ref() {
            for i in 0..bernoulli.len() {
                if let Some(breakpoint) = bernoulli.generate(i, rng) {
                    self.breakpoints.push(breakpoint);
                }
            }
        }
        let uint = rand_distr::Uniform::new(0., 1.);
        for i in self.independent_assortment.iter() {
            if rng.sample(uint) < 0.5 {
                self.breakpoints.push(i.at());
            } 
        }
        self.breakpoints.sort();
    }
    fn breakpoints(&self) -> &[Breakpoint] {
        &self.breakpoints
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MyGeneticMap {}

    impl GenerateBreakpoints for MyGeneticMap {
        fn generate_breakpoints<T: Rng>(&mut self, _rng: &mut T) {}
        fn breakpoints(&self) -> &[Breakpoint] {
            &[]
        }
    }
}

#[test]
fn test_poisson_crossover() {
    assert!(PoissonCrossover::new(0, 1, 1e-3).is_some());
    assert!(PoissonCrossover::new(0, 1, -1e-3).is_none());
    assert!(PoissonCrossover::new(0, 1, f64::NEG_INFINITY).is_none());
    assert!(PoissonCrossover::new(0, 1, f64::NAN).is_none());
    assert!(PoissonCrossover::new(1, 0, 1e-3).is_none());
}

#[test]
fn test_binomial_crossover() {
    assert!(BernoulliCrossover::new(0, 1, 1e-3).is_some());
    assert!(BernoulliCrossover::new(0, 1, 1.0).is_some());
    assert!(BernoulliCrossover::new(0, 1, 1.0 + f64::EPSILON).is_none());
    assert!(BernoulliCrossover::new(0, 1, -1e-3).is_none());
    assert!(BernoulliCrossover::new(0, 1, f64::NEG_INFINITY).is_none());
    assert!(BernoulliCrossover::new(0, 1, f64::NAN).is_none());
    assert!(BernoulliCrossover::new(1, 0, 1e-3).is_none());
}

#[test]
fn test_independent_assortment() {
    assert!(IndependentAssortment::new(1).is_some());
    assert!(IndependentAssortment::new(0).is_some());
    assert!(IndependentAssortment::new(-1).is_none());
}

#[test]
fn test_builder() {
    let builder = GeneticMapBuilder::default()
        .extend_poisson(&[PoissonCrossover::new(0, 1, 1e-3).unwrap()])
        .extend_bernoulli(&[BernoulliCrossover::new(0, 1, 0.25).unwrap()])
        .extend_poisson(&[PoissonCrossover::new(1, 2, 2e-3).unwrap()])
        .extend_independent_assortment(&[IndependentAssortment::new(1).unwrap()]);
    assert_eq!(builder.poisson().len(), 2);
    assert_eq!(builder.bernoulli().len(), 1);
    assert_eq!(builder.independent_assortment().len(), 1);
}

#[test]
fn test_breakpoint() {
    assert!(Breakpoint::Crossover(10.into()) < Breakpoint::IndependentAssortment(11.into()));
    assert!(Breakpoint::Crossover(10.into()) < Breakpoint::Crossover(11.into()));
    assert!(
        Breakpoint::IndependentAssortment(10.into()) < Breakpoint::IndependentAssortment(11.into())
    );
    assert!(Breakpoint::IndependentAssortment(10.into()) < Breakpoint::Crossover(11.into()));
}
