pub trait RowId: Sized {
    // The low-level type for the trait
    type LLType;
}
