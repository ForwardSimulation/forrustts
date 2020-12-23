use crate::nested_forward_list::NestedForwardList;
use crate::segment::Segment;

/// Data type used for edge buffering.
/// Simplification of simulated data happens
/// via [``crate::simplify_from_edge_buffer()``].
///
/// # Example
///
/// For a full example of use in simulation,
/// see the source code for
/// [``crate::wright_fisher::neutral_wf()``].
pub type EdgeBuffer = NestedForwardList<Segment>;
