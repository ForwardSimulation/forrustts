//! Core types for genetics

mod genetic_maps;

pub use genetic_maps::BernoulliCrossover;
pub use genetic_maps::Breakpoint;
pub use genetic_maps::GenerateBreakpoints;
pub use genetic_maps::GeneticMap;
pub use genetic_maps::GeneticMapBuilder;
pub use genetic_maps::GeneticMapStatus;
pub use genetic_maps::IndependentAssortment;
pub use genetic_maps::PoissonCrossover;

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

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
    assert_eq!(builder.validate(), GeneticMapStatus::Valid);
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

#[test]
fn test_breakpoint_sorting() {
    let mut v = vec![
        Breakpoint::Crossover(20.into()),
        Breakpoint::IndependentAssortment(30.into()),
        Breakpoint::IndependentAssortment(11.into()),
        Breakpoint::Crossover(10.into()),
    ];
    v.sort();
    assert_eq!(
        &v,
        &[
            Breakpoint::Crossover(10.into()),
            Breakpoint::IndependentAssortment(11.into()),
            Breakpoint::Crossover(20.into()),
            Breakpoint::IndependentAssortment(30.into()),
        ]
    );
}
