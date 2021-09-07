#![feature(test)]
extern crate test;

use forrustts_rng::Rng;
use forrustts_tables_trees::Position;

pub trait GeneticMapElement {
    //fn begin(&self) -> Position;
    //fn end(&self) -> Position;
    // Should we expect a &[Position] instead?
    fn generate_breakpoints(&self, rng: &mut Rng, breakpoints: &mut Vec<Position>);
}

pub trait GeneticMapElement2 {
    //fn begin(&self) -> Position;
    //fn end(&self) -> Position;
    // Should we expect a &[Position] instead?
    fn generate_breakpoints2(&mut self, rng: &mut Rng) -> &[Position];
}

pub struct PoissonInterval {
    beg: Position,
    end: Position,
    mean: f64,
    breakpoints: Vec<Position>,
}

impl PoissonInterval {
    pub fn new<P: Into<Position>>(beg: P, end: P, mean: f64) -> Self {
        Self {
            beg: beg.into(),
            end: end.into(),
            mean,
            breakpoints: vec![],
        }
    }
}

impl GeneticMapElement for PoissonInterval {
    fn generate_breakpoints(&self, rng: &mut Rng, breakpoints: &mut Vec<Position>) {
        use forrustts_rng::poisson;

        let n = poisson(rng, self.mean);
        for _ in 0..n {
            breakpoints.push(self.beg);
        }
    }
}

impl GeneticMapElement2 for PoissonInterval {
    fn generate_breakpoints2(&mut self, rng: &mut Rng) -> &[Position] {
        use forrustts_rng::poisson;

        self.breakpoints.clear();

        let n = poisson(rng, self.mean);
        for _ in 0..n {
            self.breakpoints.push(self.beg);
        }
        &self.breakpoints
    }
}
#[bench]
fn test2(bench: &mut test::Bencher) {
    let mut v = vec![PoissonInterval::new(0, 1, 1.25)];
    v.push(PoissonInterval::new(0, 1, 25.));
    let mut rng = Rng::new(43);
    bench.iter(|| {
        let mut b = vec![];
        for i in &mut v {
            b.extend_from_slice(i.generate_breakpoints2(&mut rng));
        }
    });
}

#[bench]
fn test1(bench: &mut test::Bencher) {
    let mut v = vec![PoissonInterval::new(0, 1, 1.25)];
    v.push(PoissonInterval::new(0, 1, 25.));
    let mut rng = Rng::new(43);
    bench.iter(|| {
        let mut b = vec![];
        for i in &mut v {
            i.generate_breakpoints(&mut rng, &mut b);
        }
    });
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
