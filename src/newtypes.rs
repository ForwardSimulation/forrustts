pub(crate) type IdType = i32;

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct NodeId(pub(crate) IdType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct EdgeId(pub(crate) IdType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct SiteId(pub(crate) IdType);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct MutationId(pub(crate) IdType);

impl_integer_ancestry_type!(NodeId, IdType, -1);
impl_integer_ancestry_type!(EdgeId, IdType, -1);
impl_integer_ancestry_type!(SiteId, IdType, -1);
impl_integer_ancestry_type!(MutationId, IdType, -1);
impl_integer_ancestry_type!(DemeId, IdType, 0);

impl_nullable_integer_ancestry_type!(NodeId);
impl_nullable_integer_ancestry_type!(EdgeId);
impl_nullable_integer_ancestry_type!(SiteId);
impl_nullable_integer_ancestry_type!(MutationId);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct DemeId(pub(crate) IdType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(pub(crate) i64);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Time(pub(crate) f64);

impl Position {
    pub const MIN: Position = Position(i64::MIN);
    pub const MAX: Position = Position(i64::MAX);

    pub fn value(&self) -> i64 {
        self.0
    }
}

impl Time {
    pub const MIN: Time = Time(f64::MIN);
    pub const MAX: Time = Time(f64::MAX);
}

impl From<i64> for Position {
    fn from(value: i64) -> Self {
        Self(value)
    }
}

impl From<f64> for Time {
    fn from(value: f64) -> Self {
        Self(value)
    }
}

impl From<i64> for Time {
    fn from(value: i64) -> Self {
        Self(value as f64)
    }
}

impl From<i32> for Time {
    fn from(value: i32) -> Self {
        Self(value as f64)
    }
}

impl PartialOrd<Time> for Time {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.0.partial_cmp(&other.0) {
            None => panic!("fatal: partial_cmp for Time received non-finite values"),
            Some(x) => Some(x),
        }
    }
}
