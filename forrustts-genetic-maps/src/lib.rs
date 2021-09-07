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

pub trait GeneticMap {
    fn generate_breakpoints(&mut self, rng: &mut Rng);
    fn breakpoints(&self) -> &[Position];
}

pub struct SimpleGeneticMap {
    map: Vec<Box<dyn GeneticMapElement>>,
    breakpoints: Vec<Position>,
}

impl SimpleGeneticMap {
    fn new() -> Self {
        let mut map: Vec<Box<dyn GeneticMapElement>> = vec![];
        map.push(Box::new(PoissonInterval::new(0, 1, 5e-2)));
        map.push(Box::new(PoissonInterval::new(0, 1, 25e-2)));
        Self {
            map,
            breakpoints: vec![],
        }
    }
}

impl GeneticMap for SimpleGeneticMap {
    fn generate_breakpoints(&mut self, rng: &mut Rng) {
        self.breakpoints.clear();
        for i in &self.map {
            i.generate_breakpoints(rng, &mut self.breakpoints);
        }
        self.breakpoints.sort();
        if !self.breakpoints.is_empty() {
            self.breakpoints.push(Position::MAX);
        }
    }

    fn breakpoints(&self) -> &[Position] {
        &self.breakpoints
    }
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
        use forrustts_rng::{poisson, uniform_i64};

        let n = poisson(rng, self.mean);
        for _ in 0..n {
            breakpoints.push(uniform_i64(rng, self.beg.into(), self.end.into()).into());
        }
    }
}

impl GeneticMapElement2 for PoissonInterval {
    fn generate_breakpoints2(&mut self, rng: &mut Rng) -> &[Position] {
        use forrustts_rng::{poisson, uniform_i64};

        self.breakpoints.clear();

        let n = poisson(rng, self.mean);
        for _ in 0..n {
            self.breakpoints
                .push(uniform_i64(rng, self.beg.into(), self.end.into()).into());
        }
        &self.breakpoints
    }
}
#[bench]
fn test2(bench: &mut test::Bencher) {
    let mut v = vec![PoissonInterval::new(0, 1, 5e-2)];
    v.push(PoissonInterval::new(0, 1, 25e-2));
    let mut rng = Rng::new(43);
    bench.iter(|| {
        for _ in 0..1000 {
            let mut b = vec![];
            for i in &mut v {
                b.extend_from_slice(i.generate_breakpoints2(&mut rng));
            }
            b.sort();
        }
    });
}

#[bench]
fn test0(bench: &mut test::Bencher) {
    let mut m = SimpleGeneticMap::new();
    let mut rng = Rng::new(43);
    bench.iter(|| {
        for _ in 0..1000 {
            m.generate_breakpoints(&mut rng);
        }
    });
}

#[bench]
fn test1(bench: &mut test::Bencher) {
    let mut v = vec![PoissonInterval::new(0, 1, 5e-2)];
    v.push(PoissonInterval::new(0, 1, 25e-2));
    let mut rng = Rng::new(43);
    bench.iter(|| {
        for _ in 0..1000 {
            let mut b = vec![];
            for i in &mut v {
                i.generate_breakpoints(&mut rng, &mut b);
            }
            b.sort();
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
