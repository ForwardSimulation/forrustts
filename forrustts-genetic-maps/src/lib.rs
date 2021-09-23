use forrustts_rng::Rng;
use forrustts_tables_trees::Position;

pub trait GeneticMapElement {
    fn begin(&self) -> Position;
    fn end(&self) -> Position;
    fn generate_breakpoints(&self, rng: &mut Rng, breakpoints: &mut Vec<Position>);
}

pub trait GeneticMap {
    fn generate_breakpoints(&mut self, rng: &mut Rng);
    fn breakpoints(&self) -> &[Position];
}

pub struct BoxedGeneticMap {
    map: Vec<Box<dyn GeneticMapElement>>,
    breakpoints: Vec<Position>,
}

impl BoxedGeneticMap {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn new_from_vec(map: Vec<Box<dyn GeneticMapElement>>) -> Self {
        Self {
            map,
            breakpoints: vec![],
        }
    }

    pub fn add_element<E: GeneticMapElement + 'static>(&mut self, element: E) {
        self.map.push(Box::new(element));
    }
}

impl Default for BoxedGeneticMap {
    fn default() -> Self {
        Self {
            map: vec![],
            breakpoints: vec![],
        }
    }
}

impl GeneticMap for BoxedGeneticMap {
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

#[derive(Copy, Clone)]
pub struct PoissonInterval {
    beg: Position,
    end: Position,
    mean: f64,
}

impl PoissonInterval {
    pub fn new<P: Into<Position>>(beg: P, end: P, mean: f64) -> Self {
        Self {
            beg: beg.into(),
            end: end.into(),
            mean,
        }
    }
}

impl GeneticMapElement for PoissonInterval {
    fn begin(&self) -> Position {
        self.beg
    }

    fn end(&self) -> Position {
        self.end
    }

    fn generate_breakpoints(&self, rng: &mut Rng, breakpoints: &mut Vec<Position>) {
        use forrustts_rng::{poisson, uniform_i64};

        let n = poisson(rng, self.mean);
        for _ in 0..n {
            breakpoints.push(uniform_i64(rng, self.beg.into(), self.end.into()).into());
        }
    }
}

#[derive(Copy, Clone)]
pub struct BinomialPoint {
    position: Position,
    probability: f64,
}

impl BinomialPoint {
    pub fn new<P: Into<Position>>(position: P, probability: f64) -> Self {
        Self {
            position: position.into(),
            probability,
        }
    }
}

impl GeneticMapElement for BinomialPoint {
    fn begin(&self) -> Position {
        self.position
    }

    fn end(&self) -> Position {
        self.position
    }

    fn generate_breakpoints(&self, rng: &mut Rng, breakpoints: &mut Vec<Position>) {
        use forrustts_rng::bernoulli_trial;

        if bernoulli_trial(rng, self.probability) {
            breakpoints.push(self.position);
        }
    }
}

#[derive(Copy, Clone)]
pub struct BinomialInterval {
    beg: Position,
    end: Position,
    probability: f64,
}

impl BinomialInterval {
    pub fn new<P: Into<Position>>(beg: P, end: P, probability: f64) -> Self {
        Self {
            beg: beg.into(),
            end: end.into(),
            probability,
        }
    }
}

impl GeneticMapElement for BinomialInterval {
    fn begin(&self) -> Position {
        self.beg
    }

    fn end(&self) -> Position {
        self.end
    }

    fn generate_breakpoints(&self, rng: &mut Rng, breakpoints: &mut Vec<Position>) {
        use forrustts_rng::{bernoulli_trial, uniform_i64};

        if bernoulli_trial(rng, self.probability) {}
        breakpoints.push(uniform_i64(rng, self.beg.into(), self.end.into()).into());
    }
}
