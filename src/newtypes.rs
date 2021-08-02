#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct NodeId(i32);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct DemeId(i32);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct EdgeId(i32);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct SiteId(i32);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct MutationId(i32);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(i64);

#[repr(transparent)]
#[derive(Copy, Clone, Debug)]
pub struct Time(f64);

/// "Null" identifier value for [``NodeId``]
pub const NULL_NODE_ID: NodeId = NodeId(-1);

impl_int_id_traits!(NodeId, i32);
impl_int_id_traits!(DemeId, i32);
impl_int_id_traits!(EdgeId, i32);
impl_int_id_traits!(SiteId, i32);
impl_int_id_traits!(MutationId, i32);
impl_int_id_traits!(Position, i64);
