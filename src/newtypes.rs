type LowLevelIdType = i32;
// TODO: remove
pub(crate) type IdType = LowLevelIdType;

pub type NodeIdLLType = LowLevelIdType;
pub type EdgeIdLLType = LowLevelIdType;
pub type MutationIdLLType = LowLevelIdType;
pub type SiteIdLLType = LowLevelIdType;
pub type DemeIdLLType = LowLevelIdType;
pub type PositonLLType = i64;
pub type TimeLLType = f64;

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct NodeId(pub(crate) NodeIdLLType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct EdgeId(pub(crate) EdgeIdLLType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct SiteId(pub(crate) SiteIdLLType);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct MutationId(pub(crate) MutationIdLLType);

impl_integer_ancestry_type!(NodeId, NodeIdLLType, -1);
impl_integer_ancestry_type!(EdgeId, EdgeIdLLType, -1);
impl_integer_ancestry_type!(SiteId, SiteIdLLType, -1);
impl_integer_ancestry_type!(MutationId, MutationIdLLType, -1);
impl_integer_ancestry_type!(DemeId, DemeIdLLType, 0);

impl_nullable_integer_ancestry_type!(NodeId);
impl_nullable_integer_ancestry_type!(EdgeId);
impl_nullable_integer_ancestry_type!(SiteId);
impl_nullable_integer_ancestry_type!(MutationId);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct DemeId(pub(crate) DemeIdLLType);

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(pub(crate) PositonLLType);

#[repr(transparent)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Time(pub(crate) TimeLLType);

impl Position {
    pub const MIN: Position = Position(PositonLLType::MIN);
    pub const MAX: Position = Position(PositonLLType::MAX);

    pub fn value(&self) -> i64 {
        self.0
    }
}

impl Time {
    pub const MIN: Time = Time(TimeLLType::MIN);
    pub const MAX: Time = Time(TimeLLType::MAX);
}

impl From<PositonLLType> for Position {
    fn from(value: PositonLLType) -> Self {
        Self(value)
    }
}

impl From<TimeLLType> for Time {
    fn from(value: TimeLLType) -> Self {
        Self(value)
    }
}

impl From<i64> for Time {
    fn from(value: i64) -> Self {
        Self(value as TimeLLType)
    }
}

impl From<i32> for Time {
    fn from(value: i32) -> Self {
        Self(value as TimeLLType)
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
