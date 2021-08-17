#![macro_use]

macro_rules! iterator_for_nodeiterator {
    ($ty: ty) => {
        impl Iterator for $ty {
            type Item = $crate::newtypes::NodeId;
            fn next(&mut self) -> Option<Self::Item> {
                self.next_node();
                self.current_node()
            }
        }
    };
}

macro_rules! impl_table_id {
    ($idtype: ident, $integer_type: ty) => {
        impl $idtype {
            pub const NULL: $idtype = Self(-1);
        }

        impl $crate::traits::private_traits::TableIdPrivate for $idtype {
            fn new(value: $crate::newtypes::TablesIdInteger) -> Self {
                Self(value)
            }

            fn new_null() -> Self {
                Self(-1)
            }

            fn raw(&self) -> $crate::newtypes::TablesIdInteger {
                self.0
            }
        }

        impl $crate::traits::private_traits::TableTypePrivate for $idtype {}

        impl $crate::traits::TableType for $idtype {}

        impl $crate::traits::TableId for $idtype {
            fn is_null(&self) -> bool {
                *self == Self::NULL
            }
        }

        impl $crate::traits::TableTypeIntoRaw for $idtype {
            type RawType = $integer_type;
            fn into_raw(self) -> Self::RawType {
                self.0
            }
        }

        impl From<$integer_type> for $idtype {
            fn from(value: $integer_type) -> Self {
                if value >= 0 {
                    Self(value)
                } else {
                    Self::NULL
                }
            }
        }

        impl From<usize> for $idtype {
            fn from(value: usize) -> Self {
                use num_traits::ToPrimitive;
                use $crate::traits::private_traits::TableIdPrivate;
                match value.to_i32() {
                    Some(x) => Self::new(x),
                    None => Self::NULL,
                }
            }
        }

        impl From<$idtype> for usize {
            fn from(value: $idtype) -> Self {
                value.0 as Self
            }
        }

        impl From<$idtype> for $crate::newtypes::TablesIdInteger {
            fn from(item: $idtype) -> Self {
                item.0
            }
        }

        impl From<i64> for $idtype {
            fn from(value: i64) -> Self {
                use num_traits::ToPrimitive;
                use $crate::traits::private_traits::TableIdPrivate;
                match value.to_i32() {
                    Some(x) => Self::new(num_traits::clamp(
                        x,
                        -1,
                        $crate::newtypes::TablesIdInteger::MAX,
                    )),
                    None => Self::NULL,
                }
            }
        }

        impl PartialEq<$integer_type> for $idtype {
            fn eq(&self, other: &$integer_type) -> bool {
                self.0 == *other
            }
        }

        impl PartialEq<$idtype> for $integer_type {
            fn eq(&self, other: &$idtype) -> bool {
                *self == other.0
            }
        }

        impl PartialOrd<$integer_type> for $idtype {
            fn partial_cmp(&self, other: &$integer_type) -> Option<std::cmp::Ordering> {
                self.0.partial_cmp(other)
            }
        }

        impl PartialOrd<$idtype> for $integer_type {
            fn partial_cmp(&self, other: &$idtype) -> Option<std::cmp::Ordering> {
                self.partial_cmp(&other.0)
            }
        }
    };
}
