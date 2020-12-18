use crate::tsdef::{IdType, Position};

#[derive(Clone, Copy)]
pub struct Segment {
    pub left: Position,
    pub right: Position,
    pub node: IdType,
}
