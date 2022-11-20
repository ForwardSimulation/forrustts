//! Compact representation of multiple forward linked lists.
//!
//! This module defines [``NestedForwardList``].  Client
//! code will typically use this type via [``crate::EdgeBuffer``].
//! See the documentation of that type for details.
//!
//! Most of API for this type is used internally, but
//! it is public in case anyone finds other uses for
//! this data structure.

use thiserror::Error;

/// Errror type for [``NestedForwardList``] operations.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum NestedForwardListError {
    /// Tail of a list is unexpectedly null.
    #[error("Tail is null")]
    NullTail,
    /// Used for invalid index values.
    #[error("Invalid index")]
    InvalidIndex,
}

/// The type used to retrieve data from [`NestedForwardList`].
type IndexType = i32;

/// Result type for [``NestedForwardList``] operations.
pub type Result<T> = std::result::Result<T, NestedForwardListError>;

/// The null value for an index
pub const NULL_INDEX: IndexType = -1;

struct ValueIterator<'list, Value> {
    list: &'list NestedForwardList<Value>,
    current: IndexType,
}

impl<'list, Value> Iterator for ValueIterator<'list, Value> {
    type Item = &'list Value;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current != NULL_INDEX {
            let rv = self.list.data_.get(self.current as usize);
            self.current = self.list.next_[self.current as usize];
            rv
        } else {
            None
        }
    }
}

struct ListIndexIterator<'list, Value> {
    list: &'list NestedForwardList<Value>,
    current: IndexType,
}

impl<'list, Value> Iterator for ListIndexIterator<'list, Value> {
    type Item = IndexType;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current != NULL_INDEX {
            let rv = self.current;
            self.current = self.list.next_[self.current as usize];
            Some(rv)
        } else {
            None
        }
    }
}

// NOTE: I am unclear how to add a Key
// to this generic.  Unlike C++, there's
// no notion of a static_assert.  Gotta Google
// this in the future!

/// Representation of multiple forward linked
/// lists flattend into vectors.
///
/// # Overview
///
/// A typical representation of a forward list
/// involves `head`, `next`, and `tail` pointers
/// to allow iteration in one direction over
/// a `Value`.
///
/// The memory layout of such a list could be
/// "node" pointers allocated on the heap.
/// For example, the node pointer could be intrusive,
/// responsible for managing ownership of the `Value`
/// elements.
///
/// Alternately, we could use flat arrays.
///
/// This type is capable of holding multiple
/// linked lists in a single set of flattened arrays.
///
/// Internally, the `head` array is accessed by list indexes.
/// In other words, to access the i-th list, make a request
/// for the i-th `head`.  Likewise for `tail`.
///
pub struct NestedForwardList<Value> {
    head_: Vec<IndexType>,
    tail_: Vec<IndexType>,
    next_: Vec<IndexType>,
    data_: Vec<Value>,
}

impl<Value> NestedForwardList<Value> {
    fn insert_new_record(&mut self, k: IndexType, v: Value) {
        self.data_.push(v);
        let x = (self.data_.len() - 1) as i32;
        self.head_[k as usize] = x;
        self.tail_[k as usize] = self.head_[k as usize];
        self.next_.push(NULL_INDEX);
    }

    fn check_key(&self, k: IndexType) -> Result<()> {
        if k < 0 {
            Err(NestedForwardListError::InvalidIndex)
        } else {
            Ok(())
        }
    }

    fn check_key_range(&self, k: usize, n: usize) -> Result<()> {
        if k >= n {
            Err(NestedForwardListError::InvalidIndex)
        } else {
            Ok(())
        }
    }

    // Public functions:

    /// Create a new instance
    pub const fn new() -> NestedForwardList<Value> {
        NestedForwardList {
            head_: Vec::<IndexType>::new(),
            tail_: Vec::<IndexType>::new(),
            next_: Vec::<IndexType>::new(),
            data_: Vec::<Value>::new(),
        }
    }

    /// Add an element to a list starting at index `k`.
    pub fn extend(&mut self, k: IndexType, v: Value) -> Result<()> {
        self.check_key(k)?;

        let idx = k as usize;
        if idx >= self.head_.len() {
            self.head_.resize(idx + 1, NULL_INDEX);
            self.tail_.resize(idx + 1, NULL_INDEX);
        }

        if self.head_[idx] == NULL_INDEX {
            self.insert_new_record(idx as IndexType, v);
            return Ok(());
        }
        let t = self.tail_[idx];
        if t == NULL_INDEX {
            return Err(NestedForwardListError::NullTail {});
        }
        self.data_.push(v);
        self.tail_[idx] = (self.data_.len() - 1) as IndexType;
        // NOTE: this differs from fwdpp, to avoid cast,
        // and thus could be failure point?
        self.next_[t as usize] = self.tail_[idx];
        self.next_.push(NULL_INDEX);
        Ok(())
    }

    // Is index at indeed part of list i?
    // Complexity: linear in length of list i
    #[allow(dead_code)]
    pub(crate) fn validate_presence_in_list(&self, i: IndexType, at: IndexType) -> bool {
        for idx in self.list_indexes_iter(i) {
            if idx == at {
                return true;
            }
        }
        false
    }

    // Insert v into list i such that it is next after after in data_
    // Corner case: what if this affects the tail?
    // NOTE: there is no check that prev and at are indeed in list i.
    #[allow(dead_code)]
    pub(crate) fn insert_after(&mut self, i: IndexType, after: IndexType, v: Value) -> Result<()> {
        self.check_key(after)?;
        self.check_key_range(after as usize, self.data_.len())?;
        debug_assert!(self.validate_presence_in_list(i, after));

        let current_next = self.next_[after as usize];

        self.data_.push(v);
        self.next_[after as usize] = (self.data_.len() as IndexType) - 1;
        self.next_.push(current_next);

        if self.tail_[i as usize] == after {
            self.tail_[i as usize] = (self.data_.len() as IndexType) - 1;
        }

        Ok(())
    }

    // For list i, insert v before before, updating next[prev]
    // if prev is NULL_INDEX, before must be head(i)
    // Complexity: amortized O(1)
    #[allow(dead_code)]
    pub(crate) fn insert_before(
        &mut self,
        i: IndexType,
        prev: IndexType,
        before: IndexType,
        v: Value,
    ) -> Result<()> {
        self.check_key(before)?;
        self.check_key_range(before as usize, self.data_.len())?;
        debug_assert!(self.validate_presence_in_list(i, before));

        self.data_.push(v);

        if prev == NULL_INDEX {
            // Then we are inserting before the head
            debug_assert_eq!(self.head_[i as usize], before);
            self.head_[i as usize] = (self.data_.len() as IndexType) - 1;
        } else {
            self.check_key(prev)?;
            self.check_key_range(prev as usize, self.data_.len())?;
            debug_assert!(self.validate_presence_in_list(i, prev));
            self.next_[prev as usize] = (self.data_.len() as IndexType) - 1;
        }
        self.next_.push(before);

        Ok(())
    }

    // at must be next[prev]
    // update data/next so that next[prev] is next[at], thus removing at
    // from the list
    // returns the new next from prev
    // Corner case: what if this affects the tail?
    // Uh, oh, we don't know how to get the tail!!!
    // The solution is clunky, requiring multiple bits of info.
    // NOTE: there is no check that prev and at are indeed in list i.
    // NOTE: the above note is not true for debug mode :)
    // An alternate solution is to have a back-index for all data elements
    // mapping them back to their parent list.
    #[allow(dead_code)]
    pub(crate) fn drop_value_from_list_at(
        &mut self,
        i: IndexType,
        prev: IndexType,
        at: IndexType,
    ) -> Result<IndexType> {
        self.check_key(prev)?;
        self.check_key(at)?;
        self.check_key_range(prev as usize, self.data_.len())?;
        self.check_key_range(at as usize, self.data_.len())?;
        debug_assert!(self.validate_presence_in_list(i, prev));
        debug_assert!(self.validate_presence_in_list(i, at));
        assert_eq!(self.next_[prev as usize], at); // TODO: convert to error
        self.next_[prev as usize] = self.next_[at as usize];
        self.next_[at as usize] = NULL_INDEX;

        // If we've dropped a tail value,
        // we need to update
        if at == self.tail_[i as usize] {
            self.tail_[i as usize] = prev;
        }

        Ok(self.next_[prev as usize])
    }

    /// Get a mutable reference to a `Value`.
    ///
    /// See [``NestedForwardList::fetch``]
    /// to get an immutable reference.
    #[inline]
    pub fn fetch_mut(&mut self, at: IndexType) -> Result<&mut Value> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.data_.len())?;
        Ok(&mut self.data_[at as usize])
    }

    /// Get the index of the head entry of a list
    /// beginning at index `at`.
    #[inline]
    pub fn head(&self, at: IndexType) -> Result<IndexType> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.head_.len())?;
        Ok(self.head_[at as usize])
    }

    /// Get the index of the tail entry
    /// for the i-th list.
    #[inline]
    pub fn tail(&self, i: IndexType) -> Result<IndexType> {
        self.check_key(i)?;
        self.check_key_range(i as usize, self.tail_.len())?;
        Ok(self.tail_[i as usize])
    }

    /// Get the index of the next data element in a list
    #[inline]
    pub fn next(&self, at: IndexType) -> Result<IndexType> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.next_.len())?;
        Ok(self.next_[at as usize])
    }

    /// Clears all data.
    /// Memory is not released.
    pub fn clear(&mut self) {
        self.data_.clear();
        self.head_.clear();
        self.tail_.clear();
        self.next_.clear();
    }

    /// Set the head/tail elements of the list
    /// beginning at index ``at`` to [``NestedForwardList::null``].
    ///
    /// # Notes
    ///
    /// This effectively "kills off" the list, preventing traversal.
    /// However, the internal data are unaffected and there is
    /// no attempt to reclaim memory.
    ///
    /// A future release may implement a stack of nullified locations,
    /// allowing their re-use.
    ///
    pub fn nullify_list(&mut self, at: IndexType) -> Result<()> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.head_.len())?;
        self.check_key_range(at as usize, self.tail_.len())?;
        self.head_[at as usize] = NULL_INDEX;
        self.tail_[at as usize] = NULL_INDEX;
        Ok(())
    }

    /// Executes the following:
    ///
    /// 1. Clears all data via [``NestedForwardList::clear``].
    /// 2. Sets the size of the number of lists to `newize`.
    /// 3. The lists are all empty, and this have head/tail
    ///    values of [``NestedForwardList::null``].
    pub fn reset(&mut self, newsize: usize) {
        self.clear();
        self.head_.resize(newsize, NULL_INDEX);
        self.tail_.resize(newsize, NULL_INDEX);
        // .fill() is experimental right now...
        //self.head_.fill(NULL_INDEX);
        //self.tail_.fill(NULL_INDEX);
        // ... so we do this lambda mapping instead
        self.head_.iter_mut().for_each(|x| *x = NULL_INDEX);
        self.tail_.iter_mut().for_each(|x| *x = NULL_INDEX);
    }

    /// Return an [`Iterator`] over the values in a given list.
    ///
    /// # Parameters
    ///
    /// * `i`: The index of a list.
    ///        An iterator is returned over the values of the `i-th` list.
    ///
    /// # Panics
    ///
    /// If `i` is out of range.
    pub fn values_iter(&self, i: IndexType) -> impl Iterator<Item = &Value> + '_ {
        ValueIterator {
            list: self,
            current: self.head(i).unwrap(),
        }
    }

    /// Return an [`Iterator`] over the indexes of all elements in the i-th list.
    ///
    /// The indexes can return the values via [`NestedForwardList::fetch`].
    /// For direct iteration over values, use [`NestedForwardList::values_iter`].
    pub fn list_indexes_iter(&self, i: IndexType) -> impl Iterator<Item = IndexType> + '_ {
        ListIndexIterator {
            list: self,
            current: self.head(i).unwrap(),
        }
    }

    /// Return an iterator over the reversed head indexes
    pub fn index_rev(&self) -> std::iter::Rev<std::ops::Range<IndexType>> {
        (0..self.len() as IndexType).rev()
    }

    /// Return length of the vector
    /// of list heads.
    pub fn len(&self) -> usize {
        self.head_.len()
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    type ListType = NestedForwardList<i32>;

    struct Datum {
        datum: i32,
    }

    type DatumList = NestedForwardList<Datum>;

    fn make_data_for_testing() -> ListType {
        let mut list = ListType::new();
        list.reset(2);
        for i in 0..3 {
            list.extend(0, 2 * i).unwrap();
        }
        for i in 0..5 {
            list.extend(1, 3 * i).unwrap();
        }
        assert_eq!(list.head_.len(), 2);
        assert_eq!(list.data_.len(), 8);
        assert_eq!(list.tail_.len(), 2);
        list
    }

    #[test]
    fn test_head_tail() {
        let mut list = make_data_for_testing();
        assert_eq!(*list.fetch_mut(list.head(0).unwrap()).unwrap(), 0);
        assert_eq!(*list.fetch_mut(list.tail(0).unwrap()).unwrap(), 4);
        assert_eq!(*list.fetch_mut(list.head(1).unwrap()).unwrap(), 0);
        assert_eq!(*list.fetch_mut(list.tail(1).unwrap()).unwrap(), 12);
    }

    #[test]
    fn test_fetch_mut() {
        let mut list = make_data_for_testing();
        let x = list.tail(1).unwrap();
        let y = list.fetch_mut(x).unwrap();
        *y += 1;
        assert_eq!(*list.fetch_mut(list.tail(1).unwrap()).unwrap(), 13);
    }

    #[test]
    fn test_fetch_mut_struct() {
        let mut list = DatumList::new();
        list.reset(1);
        list.extend(0, Datum { datum: 0 }).unwrap();
        let x = list.tail(0).unwrap();
        let y = list.fetch_mut(x).unwrap();
        y.datum = 111;
        assert_eq!(list.fetch_mut(list.tail(0).unwrap()).unwrap().datum, 111);
    }

    #[test]
    fn test_explicit_data_round_trip() {
        let mut list = make_data_for_testing();

        let mut output = Vec::<i32>::new();
        let mut itr = list.head(0).unwrap();
        while itr != NULL_INDEX {
            let val = list.fetch_mut(itr).unwrap();
            output.push(*val);
            itr = list.next(itr).unwrap();
        }
        assert_eq!(3, output.len());

        for (idx, val) in output.iter().enumerate().take(3) {
            assert_eq!(2 * idx, *val as usize);
        }

        output.clear();

        itr = list.head(1).unwrap();
        while itr != NULL_INDEX {
            let val = list.fetch_mut(itr).unwrap();
            output.push(*val);
            itr = list.next(itr).unwrap();
        }
        assert_eq!(5, output.len());

        for (idx, val) in output.iter().enumerate().take(5) {
            assert_eq!(3 * idx, *val as usize);
        }
    }

    #[test]
    fn test_nullify() {
        let mut list = make_data_for_testing();
        list.nullify_list(0).unwrap();
        let mut output = Vec::<i32>::new();
        for x in list.values_iter(0) {
            output.push(*x);
        }
        assert_eq!(output.len(), 0);
    }

    #[test]
    fn test_check_key() {
        let mut list = make_data_for_testing();
        list.reset(1);
        let result = list.extend(-1, 2);
        match result {
            Ok(_) => panic!(),
            Err(NestedForwardListError::InvalidIndex) => (),
            Err(NestedForwardListError::NullTail) => panic!(),
        }
    }

    #[test]
    fn test_head_reverse_iteration() {
        let list = make_data_for_testing();
        let mut output = Vec::<i32>::new();
        for i in list.index_rev() {
            for x in list.values_iter(i) {
                output.push(*x);
            }
        }
        assert_eq!(output.len(), 8);
        for (idx, val) in output.iter().enumerate().take(5) {
            assert_eq!(3 * idx, *val as usize);
        }
        for (idx, val) in output.iter().enumerate().skip(5) {
            assert_eq!(2 * (idx - 5), *val as usize);
        }
    }

    #[test]
    fn test_insert_after() {
        let mut list = make_data_for_testing();
        let h = list.head(1).unwrap();
        list.insert_after(1, h, 77).unwrap();
        for val in list.values_iter(1).skip(1).take(1) {
            assert_eq!(*val, 77);
        }
        let t = list.tail(1).unwrap();
        list.insert_after(1, t, -77).unwrap();
        let mut v = -1;
        for val in list.values_iter(1) {
            v = *val;
        }
        assert_eq!(v, -77);
        let t = list.tail(1).unwrap();
        assert_eq!(*list.fetch_mut(t).unwrap(), -77);
    }
}
