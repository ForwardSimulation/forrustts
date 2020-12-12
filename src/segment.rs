use crate::tsdef::TsInt;

#[derive(Clone, Copy)]
pub struct Segment {
    pub left: i64,
    pub right: i64,
    pub node: TsInt,
}
