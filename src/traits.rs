pub trait RowId: Sized + std::fmt::Debug {
    // The low-level type for the trait
    type LLType: std::fmt::Debug;
}
