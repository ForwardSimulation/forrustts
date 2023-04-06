/// A position/coordinate within a genome
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, std::hash::Hash)]
#[repr(transparent)]
pub struct Position(i64);

impl Position {
    pub fn new(position: i64) -> Option<Self> {
        if position >= 0 {
            Some(Self(position))
        } else {
            None
        }
    }

    pub fn new_valid(position: i64) -> Self {
        Self::new(position).unwrap()
    }
}

impl PartialEq<i64> for Position {
    fn eq(&self, other: &i64) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Position> for i64 {
    fn eq(&self, other: &Position) -> bool {
        *self == other.0
    }
}

impl PartialOrd<i64> for Position {
    fn partial_cmp(&self, other: &i64) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(other)
    }
}

impl PartialOrd<Position> for i64 {
    fn partial_cmp(&self, other: &Position) -> Option<std::cmp::Ordering> {
        self.partial_cmp(&other.0)
    }
}

impl TryFrom<i64> for Position {
    type Error = crate::Error;
    fn try_from(value: i64) -> Result<Self, Self::Error> {
        Self::new(value).ok_or(crate::Error::PositionError(value))
    }
}

impl From<Position> for i64 {
    fn from(value: Position) -> Self {
        value.0
    }
}
