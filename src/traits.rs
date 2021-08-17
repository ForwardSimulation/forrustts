pub(crate) mod private_traits {
    pub trait TableIdPrivate {
        fn new(value: crate::newtypes::TablesIdInteger) -> Self;
        fn new_null() -> Self;
        fn raw(&self) -> crate::newtypes::TablesIdInteger;
    }
    pub trait TableTypePrivate {}
}

/// Marker trait for low-level types in tables.
///
/// This trait cannot be implemented for types not
/// defined in this crate:
///
/// ```compile_fail
/// impl forrustts::TableType for i32 {}
/// ```
///
/// ```compile_fail
/// struct X(f32);
/// impl forrustts::TableType for X {}
/// ```
pub trait TableType: private_traits::TableTypePrivate {}

/// An integer-like object referring to a table row.
/// Trait objects can be `NULL`, indicating that
/// there is no row associated with the object.
///
/// This trait cannot be implemented for types not
/// defined in this crate because it requires [``TableType``],
/// which cannot be defined.
pub trait TableId: std::fmt::Debug + TableType + private_traits::TableIdPrivate {
    fn is_null(&self) -> bool;
}

/// Consume a [``TableId``] and return its raw
/// underlying value.
///
/// This trait cannot be implemented for types not
/// defined in this crate because it requires [``TableType``],
/// which cannot be defined.
pub trait TableTypeIntoRaw: TableType {
    /// The raw type. For example, `i32`.
    type RawType;

    /// Return the raw value.
    fn into_raw(self) -> Self::RawType;
}
