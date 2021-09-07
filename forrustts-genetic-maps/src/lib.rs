use forrustts_rng::Rng;
use forrustts_tables_trees::Position;

pub trait GeneticMapElement {
    fn begin(&self) -> Position;
    fn end(&self) -> Position;
    // Should we expect a &[Position] instead?
    fn generate_breakpoints(&self, rng: &mut Rng, breakpoints: &mut Vec<Position>);
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
