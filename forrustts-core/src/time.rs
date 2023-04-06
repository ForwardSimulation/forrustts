/// A time value
#[repr(transparent)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, std::hash::Hash, PartialOrd, Ord)]
pub struct Time(i64);

impl From<i64> for Time {
    fn from(value: i64) -> Self {
        Self(value)
    }
}

impl From<i32> for Time {
    fn from(value: i32) -> Self {
        Self(i64::from(value))
    }
}
