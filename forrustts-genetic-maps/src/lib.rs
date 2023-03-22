//! Core types for genetic maps

use forrustts_core::newtypes::Position;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum Breakpoint {
    Crossover(Position),
    IndependentAssortment(Position),
}

pub trait GenerateBreakpoints<T> {
    fn generate_breakpoints(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Breakpoint];
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MyGeneticMap {}

    impl GenerateBreakpoints<()> for MyGeneticMap {
        fn generate_breakpoints(&mut self, _rng: &mut ()) {}
        fn breakpoints(&self) -> &[Breakpoint] {
            &[]
        }
    }
}
