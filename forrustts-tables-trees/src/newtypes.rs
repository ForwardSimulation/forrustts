/// The low-level representation
/// of a [``TableId``](crate::traits::TableId)
pub type TablesIdInteger = i32;

/// The low-level representation
/// of a [``Position``].
pub type PositionLLType = i64;

/// The low-level representation
/// of a [``Time``].
pub type TimeLLType = f64;

/// A [``TableId``](crate::traits::TableId) for a node.
#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct NodeId(pub(crate) TablesIdInteger);

/// A [``TableId``](crate::traits::TableId) for an edge.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct EdgeId(pub(crate) TablesIdInteger);

/// A [``TableId``](crate::traits::TableId) for a site.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct SiteId(pub(crate) TablesIdInteger);

/// A [``TableId``](crate::traits::TableId) for a mutation.
#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct MutationId(pub(crate) TablesIdInteger);

/// A [``TableId``](crate::traits::TableId) for a deme.
#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct DemeId(pub(crate) TablesIdInteger);

impl_table_id!(NodeId, TablesIdInteger);
impl_table_id!(EdgeId, TablesIdInteger);
impl_table_id!(SiteId, TablesIdInteger);
impl_table_id!(MutationId, TablesIdInteger);
impl_table_id!(DemeId, TablesIdInteger);

/// A position/coordinate within a genome
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(pub(crate) PositionLLType);

/// A time value
#[repr(transparent)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Time(pub(crate) TimeLLType);

impl Position {
    /// Minimum value
    pub const MIN: Position = Position(PositionLLType::MIN);
    /// Maximum value
    pub const MAX: Position = Position(PositionLLType::MAX);
}

impl Time {
    /// Minimum value
    pub const MIN: Time = Time(TimeLLType::MIN);
    /// Maximum value
    pub const MAX: Time = Time(TimeLLType::MAX);
}

impl crate::traits::private_traits::TableTypePrivate for Position {}
impl crate::traits::private_traits::TableTypePrivate for Time {}
impl crate::traits::TableType for Position {}
impl crate::traits::TableType for Time {}

impl crate::traits::TableTypeIntoRaw for Position {
    type RawType = PositionLLType;
    fn into_raw(self) -> Self::RawType {
        self.0
    }
}

impl crate::traits::TableTypeIntoRaw for Time {
    type RawType = TimeLLType;
    fn into_raw(self) -> Self::RawType {
        self.0
    }
}

impl PartialEq<PositionLLType> for Position {
    fn eq(&self, other: &PositionLLType) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Position> for PositionLLType {
    fn eq(&self, other: &Position) -> bool {
        *self == other.0
    }
}

impl PartialOrd<PositionLLType> for Position {
    fn partial_cmp(&self, other: &PositionLLType) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(other)
    }
}

impl PartialOrd<Position> for PositionLLType {
    fn partial_cmp(&self, other: &Position) -> Option<std::cmp::Ordering> {
        self.partial_cmp(&other.0)
    }
}

impl From<PositionLLType> for Position {
    fn from(value: PositionLLType) -> Self {
        Self(value)
    }
}

impl From<Position> for PositionLLType {
    fn from(value: Position) -> Self {
        value.0
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

impl From<Time> for f64 {
    fn from(value: Time) -> Self {
        value.0
    }
}

impl From<Time> for i64 {
    fn from(value: Time) -> Self {
        value.0 as Self
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
