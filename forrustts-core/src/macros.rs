#![macro_use]

macro_rules! impl_table_id {
    ($idtype: ident, $integer_type: ty) => {
        impl $idtype {
            /// NULL value for the type
            pub const NULL: $idtype = Self(-1);

            fn new(value: $integer_type) -> Self {
                Self(value)
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

        impl std::fmt::Display for $idtype {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
                write!(f, "{}", self.0)
            }
        }

        impl TryFrom<$idtype> for usize {
            type Error = $crate::Error;

            fn try_from(value: $idtype) -> Result<Self, Self::Error> {
                usize::try_from(value.0).map_err(|_| {
                    $crate::Error::ConversionError(format!("could not convert {} to usize", value))
                })
            }
        }

        impl TryFrom<&$idtype> for usize {
            type Error = $crate::Error;

            fn try_from(value: &$idtype) -> Result<Self, Self::Error> {
                usize::try_from(value.0).map_err(|_| {
                    $crate::Error::ConversionError(format!("could not convert {} to usize", value))
                })
            }
        }

        impl TryFrom<usize> for $idtype {
            type Error = $crate::Error;

            fn try_from(value: usize) -> Result<Self, Self::Error> {
                let ll = <$integer_type>::try_from(value).map_err(|_| {
                    $crate::Error::ConversionError(format!("could not convert {} to usize", value))
                })?;
                Ok($idtype::new(ll))
            }
        }

        impl From<$idtype> for $integer_type {
            fn from(item: $idtype) -> Self {
                item.0
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

macro_rules! impl_get_raw {
    ($newtype: ident, $lltype: ty) => {
        impl $newtype {
            pub fn raw(self) -> $lltype {
                self.0
            }
        }
    };
}
