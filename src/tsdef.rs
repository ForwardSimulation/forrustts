/// Integer type used to refer to a genomic coordinate/position
pub type Position = i64;
/// Time is recorded as a discrete, signed integer.
pub type Time = i64;
/// Integer type used to refer to [Node](struct.Node.html) objects.
pub type IdType = i32;
pub type IdVec = Vec<IdType>;

/// Equals -1 (minus one).
/// Primary use is to indicate a null [Node](struct.Node.html) id.
pub static NULL_ID: IdType = -1;
