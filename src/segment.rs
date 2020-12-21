use crate::tsdef::{IdType, Position};

#[derive(Clone, Copy)]
pub struct Segment {
    pub left: Position,
    pub right: Position,
    pub node: IdType,
}

impl Segment {
    pub fn new(left: Position, right: Position, node: IdType) -> Self {
        Segment { left, right, node }
    }
}
