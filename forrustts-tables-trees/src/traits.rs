pub(crate) mod private_traits {
    pub trait LowLevelTableType {
        type Type;
        fn new(value: Self::Type) -> Self;
        fn raw(&self) -> Self::Type;
    }

    pub trait NullableLowLevelTableType: LowLevelTableType {
        fn new_null() -> Self;
    }
}

/// Marker trait for low-level types in tables.
///
/// This trait cannot be implemented for types not
/// defined in this crate:
///
/// ```compile_fail
/// impl forrustts::TableType for i32 {
///     type LowLevelType = i32;
///
///     fn into_raw(self) -> Self::LowLevelType {
///         self
///     }
/// }
/// ```
///
/// ```compile_fail
/// struct X(f32);
/// impl forrustts::TableType for X {
///     type LowLevelType = f64;
///
///     fn into_raw(self) -> Self::LowLevelType {
///         self.0
///     }
/// }
/// ```
pub trait TableType: private_traits::LowLevelTableType {
    /// The underlying type
    type LowLevelType;

    /// Return the underlying value
    fn into_raw(self) -> Self::LowLevelType;
}

/// An integer-like object referring to a table row.
/// Trait objects can be `NULL`, indicating that
/// there is no row associated with the object.
///
/// This trait cannot be implemented for types not
/// defined in this crate because it requires [``TableType``],
/// which cannot be defined.
pub trait TableId: std::fmt::Debug + TableType + private_traits::NullableLowLevelTableType {
    /// Return true if `self` is equal to the
    /// type's `NULL` value.
    fn is_null(&self) -> bool;
}
