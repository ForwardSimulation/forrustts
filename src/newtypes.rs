pub(crate) type IdType = i32;

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct NodeId(pub(crate) IdType);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct DemeId(pub(crate) IdType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct EdgeId(pub(crate) IdType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct SiteId(pub(crate) IdType);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct MutationId(pub(crate) IdType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(pub(crate) i64);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Time(pub(crate) f64);

impl_row_id_traits!(NodeId, IdType);
impl_row_id_traits!(EdgeId, IdType);
impl_row_id_traits!(SiteId, IdType);
impl_row_id_traits!(MutationId, IdType);

impl PartialOrd<Time> for Time {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.0.partial_cmp(&other.0) {
            None => panic!("fatal: partial_cmp for Time received non-finite values"),
            Some(x) => Some(x),
        }
    }
}
