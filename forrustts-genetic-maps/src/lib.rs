//! Core types for genetic maps

use forrustts_core::newtypes::Position;

#[non_exhaustive]
pub enum Breakpoint {
    Crossover(Position),
    IndependentAssortment(Position),
}

pub trait GeneticMap<T> {
    fn generate_breakpoints(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Breakpoint];
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MyGeneticMap {}

    impl GeneticMap<()> for MyGeneticMap {
        fn generate_breakpoints(&mut self, _rng: &mut ()) {}
        fn breakpoints(&self) -> &[Breakpoint] {
            &[]
        }
    }
}
