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
type LowLevelTimeType = f64;
type LowLevelPositonType = i64;

/// A [``TableId``](crate::traits::TableId) for a node.
///
/// ```
/// use forrustts_tables_trees::prelude::*;
///
/// let n = NodeId::from(-1);
/// assert_eq!(n, -1);
/// let r = n.into_raw();
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

pub struct TableIdRange<T>
where
    T: crate::TableId + Copy + std::cmp::Ord + std::ops::Add<Output = T> + std::ops::AddAssign,
{
    current: T,
    last: T,
    increment: T,
}

impl<T> TableIdRange<T>
where
    T: crate::TableId + Copy + std::cmp::Ord + std::ops::Add<Output = T> + std::ops::AddAssign,
{
    pub fn new(from: T, to: T, inclusive: bool) -> Self {
        let increment = if from < to { T::new(1) } else { T::new(-1) };
        let last = if inclusive { to + increment } else { to };
        Self {
            current: from,
            last,
            increment,
        }
    }
}

pub fn make_id_range<T>(from: T, to: T, inclusive: bool) -> impl Iterator<Item=T>
where
    T: crate::TableId + Copy + std::cmp::Ord + std::ops::Add<Output = T> + std::ops::AddAssign,
{
    TableIdRange::<T>::new(from, to, inclusive)
}

impl<T> Iterator for TableIdRange<T>
where
    T: crate::TableId + Copy + std::cmp::Ord + std::ops::Add<Output = T> + std::ops::AddAssign,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current != self.last {
            let rv = self.current;
            self.current += self.increment;
            Some(rv)
        } else {
            None
        }
    }
}

/// A position/coordinate within a genome
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(pub(crate) LowLevelPositonType);

/// A time value
#[repr(transparent)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Time(pub(crate) LowLevelTimeType);

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

impl_low_level_table_type!(Time, LowLevelTimeType);
impl_low_level_table_type!(Position, LowLevelPositonType);

impl std::ops::Add for Time {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl std::ops::Sub for Time {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl std::fmt::Display for Time {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Time({})", self.0)
    }
}

impl std::ops::Add for Position {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl std::ops::Sub for Position {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Position({})", self.0)
    }
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

impl From<i64> for Time {
    fn from(value: i64) -> Self {
        Self(value as LowLevelTimeType)
    }
}

impl From<i32> for Time {
    fn from(value: i32) -> Self {
        Self(value as LowLevelTimeType)
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

#[cfg(test)]
mod test_newtypes {
    use super::*;

    #[test]
    fn test_id_range() {
        {
            let r = TableIdRange::<NodeId>::new(NodeId::from(1), NodeId::from(3), false);

            let vr: Vec<NodeId> = r.collect();

            assert!(vr == vec![NodeId::from(1), NodeId::from(2)]);
        }

        {
            let r = make_id_range(NodeId::from(1), NodeId::from(3), false);

            let vr: Vec<NodeId> = r.collect();

            assert!(vr == vec![NodeId::from(1), NodeId::from(2)]);
        }

        {
            let r = TableIdRange::<NodeId>::new(NodeId::from(1), NodeId::from(3), true);

            let vr: Vec<NodeId> = r.collect();

            assert!(vr == vec![NodeId::from(1), NodeId::from(2), NodeId::from(3)]);
        }

        {
            let r = TableIdRange::<NodeId>::new(NodeId::from(3), NodeId::from(1), false);

            let vr: Vec<NodeId> = r.collect();

            assert!(vr == vec![NodeId::from(3), NodeId::from(2)]);
        }

        {
            let r = TableIdRange::<NodeId>::new(NodeId::from(3), NodeId::from(1), true);

            let vr: Vec<NodeId> = r.collect();

            assert!(vr == vec![NodeId::from(3), NodeId::from(2), NodeId::from(1)]);
        }

        {
            let r = TableIdRange::<NodeId>::new(NodeId::from(3), NodeId::from(3), false);

            let vr: Vec<NodeId> = r.collect();

            assert!(vr == Vec::<NodeId>::default());
        }

        {
            let r = TableIdRange::<NodeId>::new(NodeId::from(3), NodeId::from(3), true);

            let vr: Vec<NodeId> = r.collect();

            assert!(vr == vec![NodeId::from(3)]);
        }
    }
}
