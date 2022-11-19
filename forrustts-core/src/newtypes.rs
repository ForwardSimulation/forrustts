/// The low-level representation
/// of a [``TableId``](crate::traits::TableId)
#[deprecated(note = "use TableType::Type instead")]
pub type TablesIdInteger = i32;

/// The low-level representation
/// of a [``Position``].
#[deprecated(note = "use TableType::Type instead")]
pub type PositionLLType = i64;

/// The low-level representation
/// of a [``Time``].
#[deprecated(note = "use TableType::Type instead")]
pub type TimeLLType = f64;

type LowLevelIdType = i32;
type LowLevelTimeType = i64;
type LowLevelPositonType = i64;

/// A [``TableId``](crate::traits::TableId) for a node.
///
/// ```
/// # use forrustts_core::newtypes::NodeId;
/// let n = NodeId::from(-1);
/// assert_eq!(n, -1);
/// let r = n.raw();
/// assert_eq!(r, -1);
/// ```
#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct NodeId(pub(crate) LowLevelIdType);

/// A [``TableId``](crate::traits::TableId) for an edge.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct EdgeId(pub(crate) LowLevelIdType);

/// A [``TableId``](crate::traits::TableId) for a site.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct SiteId(pub(crate) LowLevelIdType);

/// A [``TableId``](crate::traits::TableId) for a mutation.
#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct MutationId(pub(crate) LowLevelIdType);

/// A [``TableId``](crate::traits::TableId) for a deme.
#[repr(transparent)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
pub struct DemeId(pub(crate) LowLevelIdType);

impl_table_id!(NodeId, LowLevelIdType);
impl_table_id!(EdgeId, LowLevelIdType);
impl_table_id!(SiteId, LowLevelIdType);
impl_table_id!(MutationId, LowLevelIdType);
impl_table_id!(DemeId, LowLevelIdType);

impl_get_raw!(NodeId, LowLevelIdType);
impl_get_raw!(EdgeId, LowLevelIdType);
impl_get_raw!(SiteId, LowLevelIdType);
impl_get_raw!(MutationId, LowLevelIdType);
impl_get_raw!(DemeId, LowLevelIdType);

/// A position/coordinate within a genome
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(pub(crate) LowLevelPositonType);

/// A time value
#[repr(transparent)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, std::hash::Hash)]
pub struct Time(pub(crate) LowLevelTimeType);

impl_get_raw!(Position, LowLevelPositonType);
impl_get_raw!(Time, LowLevelTimeType);

impl Position {
    /// Minimum value
    pub const MIN: Position = Position(LowLevelPositonType::MIN);
    /// Maximum value
    pub const MAX: Position = Position(LowLevelPositonType::MAX);
}

impl Time {
    /// Minimum value
    pub const MIN: Time = Time(LowLevelTimeType::MIN);
    /// Maximum value
    pub const MAX: Time = Time(LowLevelTimeType::MAX);
}

impl PartialEq<LowLevelPositonType> for Position {
    fn eq(&self, other: &LowLevelPositonType) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Position> for LowLevelPositonType {
    fn eq(&self, other: &Position) -> bool {
        *self == other.0
    }
}

impl PartialOrd<LowLevelPositonType> for Position {
    fn partial_cmp(&self, other: &LowLevelPositonType) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(other)
    }
}

impl PartialOrd<Position> for LowLevelPositonType {
    fn partial_cmp(&self, other: &Position) -> Option<std::cmp::Ordering> {
        self.partial_cmp(&other.0)
    }
}

impl From<LowLevelPositonType> for Position {
    fn from(value: LowLevelPositonType) -> Self {
        Self(value)
    }
}

impl From<Position> for LowLevelPositonType {
    fn from(value: Position) -> Self {
        value.0
    }
}

impl From<LowLevelTimeType> for Time {
    fn from(value: LowLevelTimeType) -> Self {
        Self(value)
    }
}

impl From<i32> for Time {
    fn from(value: i32) -> Self {
        Self(value as LowLevelTimeType)
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
