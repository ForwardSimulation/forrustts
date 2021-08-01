#![macro_use]

macro_rules! iterator_for_nodeiterator {
    ($ty: ty) => {
        impl Iterator for $ty {
            type Item = $crate::IdType;
            fn next(&mut self) -> Option<Self::Item> {
                self.next_node();
                self.current_node()
            }
        }
    };
}

macro_rules! impl_int_id_traits {
    ($idtype: ty, $integer_type: ty) => {
        impl $idtype {
            pub fn is_null(&self) -> bool {
                self.0 == -1
            }
        }

        impl From<$integer_type> for $idtype {
            fn from(value: $integer_type) -> Self {
                Self(value)
            }
        }

        impl From<$idtype> for usize {
            fn from(value: $idtype) -> Self {
                value.0 as usize
            }
        }

        impl From<$idtype> for $integer_type {
            fn from(value: $idtype) -> Self {
                value.0
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
