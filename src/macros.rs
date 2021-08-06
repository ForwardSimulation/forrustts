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

macro_rules! impl_integer_ancestry_type {
    ($idtype: ident, $integer_type: ty, $minval: expr) => {
        impl $idtype {
            pub fn new(value: $integer_type) -> Result<Self, $crate::error::RowIdError<$idtype>> {
                if value < $minval {
                    Err($crate::error::RowIdError::<$idtype>::InvalidValue { value })
                } else {
                    Ok(Self(value))
                }
            }
        }

        impl $crate::traits::AncestryType for $idtype {
            type LLType = $integer_type;

            fn value(&self) -> Self::LLType {
                self.0
            }
        }

        impl From<$integer_type> for $idtype {
            fn from(value: $integer_type) -> Self {
                debug_assert!(value >= $minval);
                Self(value)
            }
        }

        impl std::convert::TryFrom<usize> for $idtype {
            type Error = $crate::error::RowIdError<$idtype>;
            fn try_from(value: usize) -> Result<Self, Self::Error> {
                if value > (<$integer_type>::MAX as usize) {
                    Err(Self::Error::ValueOverflow {
                        value: value.to_string(),
                    })
                } else {
                    Self::new(value as $integer_type)
                }
            }
        }

        impl std::convert::TryFrom<$idtype> for usize {
            type Error = $crate::error::RowIdError<$idtype>;
            fn try_from(value: $idtype) -> Result<Self, Self::Error> {
                if value < 0 {
                    Err(Self::Error::InvalidValue { value: value.0 })
                } else {
                    Ok(value.0 as usize)
                }
            }
        }

        impl std::convert::TryFrom<i64> for $idtype {
            type Error = $crate::error::RowIdError<$idtype>;
            fn try_from(value: i64) -> Result<Self, Self::Error> {
                if value > (<$integer_type>::MAX as i64) {
                    Err(Self::Error::ValueOverflow {
                        value: value.to_string(),
                    })
                } else {
                    Self::new(value as $integer_type)
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

macro_rules! impl_nullable_integer_ancestry_type {
    ($idtype: ident) => {
        impl $idtype {
            pub const NULL: $idtype = Self(-1);
        }
    };
}
