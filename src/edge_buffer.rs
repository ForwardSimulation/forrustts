use crate::segment::Segment;
use crate::nested_forward_list::NestedForwardList;

pub type EdgeBuffer = NestedForwardList<Segment>;
