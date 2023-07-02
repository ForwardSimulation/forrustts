use crate::rand::Rng;
use crate::rand_distr;

use forrustts_core::Position;

/// Breakpoint positions from crossover events
///
/// # Notes
///
/// * Comparision operations are based
///   on the stored
///   [`Position`](forrustts_core::Position).
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

    /// Number of generated breakpoints.
    fn len(&self) -> usize {
        self.breakpoints().len()
    }

    /// Checks if generated breakpoints are empty.
    fn is_empty(&self) -> bool {
        self.breakpoints().is_empty()
    }
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
        L: TryInto<Position>,
        R: TryInto<Position>,
    {
        let left = left.try_into().ok()?;
        let right = right.try_into().ok()?;
        if left < 0 || right < 0 || right <= left || !mean.is_finite() || mean < 0.0 {
            None
        } else {
            Some(Self { left, right, mean })
        }
    }

    pub fn new_point<P>(at: P, mean: f64) -> Option<Self>
    where
        P: TryInto<Position>,
    {
        let left = at.try_into().ok()?;
        let right = i64::from(left) + 1;
        Self::new(left, right, mean)
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
        L: TryInto<Position>,
        R: TryInto<Position>,
    {
        let left = left.try_into().ok()?;
        let right = right.try_into().ok()?;
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
    pub fn new_point<P>(at: P, probability: f64) -> Option<Self>
    where
        P: TryInto<Position>,
    {
        let left = at.try_into().ok()?;
        let right = i64::from(left) + 1;
        Self::new(left, right, probability)
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
    pub fn new<P: TryInto<Position>>(at: P) -> Option<Self> {
        let at = at.try_into().ok()?;
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

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum GeneticMapStatus {
    Valid,
    IndependentAssortmentWithinRegions,
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

    fn validate_independent_assortment(&self) -> bool {
        for i in &self.independent_assortment {
            for p in &self.poisson {
                if i.at > p.left && i.at < p.right {
                    return false;
                }
            }
            for b in &self.bernoulli {
                if i.at > b.left && i.at < b.right {
                    return false;
                }
            }
        }

        true
    }

    pub fn validate(&self) -> GeneticMapStatus {
        if !self.validate_independent_assortment() {
            return GeneticMapStatus::IndependentAssortmentWithinRegions;
        }
        GeneticMapStatus::Valid
    }
}

#[derive(Debug)]
struct PoissonRegions {
    regions: Vec<rand_distr::Uniform<Position>>,
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
                let u = rand_distr::Uniform::new(p.left(), p.right());
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
        Breakpoint::Crossover(rng.sample(self.regions[i]))
    }
}

#[derive(Debug)]
struct BernoulliRegions {
    regions: Vec<rand_distr::Uniform<Position>>,
    probabilities: Vec<rand_distr::Bernoulli>,
}

impl BernoulliRegions {
    fn new(binomial: Vec<BernoulliCrossover>) -> Option<Self> {
        let mut regions = vec![];
        let mut probabilities = vec![];

        for b in binomial {
            let prob = b.probability();
            if prob > 0.0 {
                let u = rand_distr::Uniform::new(b.left(), b.right());
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
            Some(Breakpoint::Crossover(rng.sample(self.regions[index])))
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
        if builder.validate() != GeneticMapStatus::Valid {
            return None;
        }
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
