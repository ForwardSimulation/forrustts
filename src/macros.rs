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

macro_rules! impl_row_id_traits {
    ($idtype: ident, $integer_type: ty) => {
        impl $idtype {
            pub fn new(value: $integer_type) -> Result<Self, $crate::ForrusttsError> {
                if value < Self::NULL {
                    Err($crate::ForrusttsError::RowIdError {
                        value: value.to_string(),
                    })
                } else {
                    Ok(Self(value))
                }
            }

            pub fn is_null(&self) -> bool {
                self.0 == -1
            }

            pub fn value(&self) -> $integer_type {
                self.0
            }

            pub const NULL: $idtype = $idtype(-1);
        }

        impl $crate::traits::RowId for $idtype {
            type LLType = $idtype;
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
