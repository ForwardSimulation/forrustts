use crate::newtypes::{NodeId, Position};

/// A segment is a half-open
/// interval of [``Position``]s
/// associated with a [``NodeId``].
///
/// This type is public primarily because
/// it is the value element of a
/// [``crate::EdgeBuffer``].
#[derive(Clone, Copy)]
pub struct Segment {
    /// Left edge of interval
    pub left: Position,
    /// Right edge of interval
    pub right: Position,
    /// The node
    pub node: NodeId,
}

impl Segment {
    /// Create a new instance.
    pub fn new(left: Position, right: Position, node: NodeId) -> Self {
        Segment { left, right, node }
    }
}
