//! Postprocessing modules for filters and genotypes.

use crate::variations::{Variant, Vars};
use std::borrow::Cow;

pub mod amplicon;
pub mod simple;
pub mod somatic;

pub(crate) enum CowVarIter<'a> {
    Borrowed(std::slice::Iter<'a, Variant>),
    Owned(std::vec::IntoIter<Variant>),
}

impl<'a> Iterator for CowVarIter<'a> {
    type Item = Cow<'a, Variant>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Borrowed(it) => it.next().map(Cow::Borrowed),
            Self::Owned(it) => it.next().map(Cow::Owned),
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Self::Borrowed(it) => it.size_hint(),
            Self::Owned(it) => it.size_hint(),
        }
    }
}

impl ExactSizeIterator for CowVarIter<'_> {
    #[inline]
    fn len(&self) -> usize {
        match self {
            Self::Borrowed(it) => it.len(),
            Self::Owned(it) => it.len(),
        }
    }
}

impl CowVarIter<'_> {
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[inline]
pub(crate) fn unpack_cow_vars<'a>(
    vars: Cow<'a, Vars>,
) -> (Option<Cow<'a, Variant>>, CowVarIter<'a>, Cow<'a, str>) {
    match vars {
        Cow::Borrowed(vars) => (
            vars.ref_variant.as_ref().map(Cow::Borrowed),
            CowVarIter::Borrowed(vars.variants.iter()),
            Cow::Borrowed(vars.sv.as_str()),
        ),
        Cow::Owned(Vars {
            ref_variant,
            variants,
            sv,
            ..
        }) => (
            ref_variant.map(Cow::Owned),
            CowVarIter::Owned(variants.into_iter()),
            Cow::Owned(sv),
        ),
    }
}
