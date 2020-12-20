use thiserror::Error;

#[derive(Error, Debug)]
pub enum NestedForwardListError {
    #[error("Invalid key")]
    InvalidKey,
    #[error("Invalid key")]
    KeyOutOfRange,
}

pub type Result<T> = std::result::Result<T, NestedForwardListError>;

// NOTE: I am unclear how to add a Key
// to this generic.  Unlike C++, there's
// no notion of a static_assert.  Gotta Google
// this in the future!

pub struct NestedForwardList<Value> {
    head_: Vec<i32>,
    tail_: Vec<i32>,
    next_: Vec<i32>,
    data_: Vec<Value>,
}

impl<Value> NestedForwardList<Value> {
    fn insert_new_record(&mut self, k: i32, v: Value) {
        self.data_.push(v);
        let x = (self.data_.len() - 1) as i32;
        self.head_[k as usize] = x;
        self.tail_[k as usize] = self.head_[k as usize];
        self.next_.push(NestedForwardList::<Value>::null());
    }

    fn check_key(&self, k: i32) -> Result<()> {
        if k < 0 {
            Err(NestedForwardListError::InvalidKey)
        } else {
            Ok(())
        }
    }

    fn check_key_range(&self, k: usize, n: usize) -> Result<()> {
        if k >= n {
            Err(NestedForwardListError::KeyOutOfRange)
        } else {
            Ok(())
        }
    }

    // Public functions:

    pub const fn new() -> NestedForwardList<Value> {
        NestedForwardList {
            head_: Vec::<i32>::new(),
            tail_: Vec::<i32>::new(),
            next_: Vec::<i32>::new(),
            data_: Vec::<Value>::new(),
        }
    }

    pub fn extend(&mut self, k: i32, v: Value) -> Result<()> {
        self.check_key(k)?;

        let idx = k as usize;
        if idx >= self.head_.len() {
            self.head_
                .resize(idx + 1, NestedForwardList::<Value>::null());
            self.tail_
                .resize(idx + 1, NestedForwardList::<Value>::null());
        }

        if self.head_[idx] == NestedForwardList::<Value>::null() {
            self.insert_new_record(idx as i32, v);
            return Ok(());
        }
        let t = self.tail_[idx];
        if t == NestedForwardList::<Value>::null() {
            panic!("unexpected null");
        }
        self.data_.push(v);
        self.tail_[idx] = (self.data_.len() - 1) as i32;
        // NOTE: this differs from fwdpp, to avoid cast,
        // and thus could be failure point?
        self.next_[t as usize] = self.tail_[idx];
        self.next_.push(NestedForwardList::<Value>::null());
        Ok(())
    }

    #[inline]
    pub fn null() -> i32 {
        -1
    }

    #[inline]
    pub fn fetch_mut(&mut self, at: i32) -> Result<&mut Value> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.data_.len())?;
        Ok(&mut self.data_[at as usize])
    }

    #[inline]
    pub fn fetch(&self, at: i32) -> Result<&Value> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.data_.len())?;
        Ok(&self.data_[at as usize])
    }

    #[inline]
    pub fn head(&self, at: i32) -> Result<i32> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.head_.len())?;
        Ok(self.head_[at as usize])
    }

    #[inline]
    pub fn tail(&self, at: i32) -> Result<i32> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.tail_.len())?;
        Ok(self.tail_[at as usize])
    }

    #[inline]
    pub fn next(&self, at: i32) -> Result<i32> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.next_.len())?;
        Ok(self.next_[at as usize])
    }

    pub fn clear(&mut self) {
        self.data_.clear();
        self.head_.clear();
        self.tail_.clear();
        self.next_.clear();
    }

    pub fn for_each(&self, at: i32, mut f: impl FnMut(&Value) -> bool) -> Result<()> {
        let mut itr = self.head(at)?;
        while itr != NestedForwardList::<Value>::null() {
            let val = self.fetch(itr)?;
            let check = f(val);
            if !check {
                break;
            }
            itr = self.next(itr)?;
        }
        Ok(())
    }

    pub fn nullify_list(&mut self, at: i32) -> Result<()> {
        self.check_key(at)?;
        self.check_key_range(at as usize, self.head_.len())?;
        self.check_key_range(at as usize, self.tail_.len())?;
        self.head_[at as usize] = NestedForwardList::<Value>::null();
        self.tail_[at as usize] = NestedForwardList::<Value>::null();
        Ok(())
    }

    pub fn reset(&mut self, newsize: usize) {
        self.clear();
        self.head_
            .resize(newsize, NestedForwardList::<Value>::null());
        self.tail_
            .resize(newsize, NestedForwardList::<Value>::null());
        // .fill() is experimental right now...
        //self.head_.fill(NestedForwardList::<Value>::null());
        //self.tail_.fill(NestedForwardList::<Value>::null());
        // ... so we do this lambda mapping instead
        self.head_.iter_mut().for_each(|x| *x = Self::null());
        self.tail_.iter_mut().for_each(|x| *x = Self::null());
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
        list
    }

    #[test]
    fn test_head_tail() {
        let list = make_data_for_testing();
        assert_eq!(*list.fetch(list.head(0).unwrap()).unwrap(), 0);
        assert_eq!(*list.fetch(list.tail(0).unwrap()).unwrap(), 4);
        assert_eq!(*list.fetch(list.head(1).unwrap()).unwrap(), 0);
        assert_eq!(*list.fetch(list.tail(1).unwrap()).unwrap(), 12);
    }

    #[test]
    fn test_fetch_mut() {
        let mut list = make_data_for_testing();
        let x = list.tail(1).unwrap();
        let y = list.fetch_mut(x).unwrap();
        *y += 1;
        assert_eq!(*list.fetch(list.tail(1).unwrap()).unwrap(), 13);
    }

    #[test]
    fn test_fetch_mut_struct() {
        let mut list = DatumList::new();
        list.reset(1);
        list.extend(0, Datum { datum: 0 }).unwrap();
        let x = list.tail(0).unwrap();
        let y = list.fetch_mut(x).unwrap();
        y.datum = 111;
        assert_eq!(list.fetch(list.tail(0).unwrap()).unwrap().datum, 111);
    }

    #[test]
    fn test_explicit_data_round_trip() {
        let list = make_data_for_testing();

        let mut output = Vec::<i32>::new();
        let mut itr = list.head(0).unwrap();
        while itr != ListType::null() {
            let val = list.fetch(itr).unwrap();
            output.push(*val);
            itr = list.next(itr).unwrap();
        }
        assert_eq!(3, output.len());

        for (idx, val) in output.iter().enumerate().take(3) {
            assert_eq!(2 * idx, *val as usize);
        }

        output.clear();

        itr = list.head(1).unwrap();
        while itr != ListType::null() {
            let val = list.fetch(itr).unwrap();
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
        list.for_each(0, |x: &i32| {
            output.push(*x);
            true
        })
        .unwrap();
        assert_eq!(output.len(), 0);
    }

    #[test]
    fn test_for_each_data_round_trip() {
        let list = make_data_for_testing();
        let mut output = Vec::<i32>::new();

        list.for_each(0, |x: &i32| {
            output.push(*x);
            true
        })
        .unwrap();

        for (idx, val) in output.iter().enumerate().take(3) {
            assert_eq!(2 * idx, *val as usize);
        }

        output.clear();

        list.for_each(1, |x: &i32| {
            output.push(*x);
            true
        })
        .unwrap();

        for (idx, val) in output.iter().enumerate().take(5) {
            assert_eq!(3 * idx, *val as usize);
        }

        output.clear();
        list.for_each(1, |_: &i32| false).unwrap();

        assert_eq!(output.len(), 0);
    }

    #[test]
    fn test_check_key() {
        let mut list = make_data_for_testing();
        list.reset(1);
        let result = list.extend(-1, 2);
        match result {
            Ok(_) => panic!(),
            Err(NestedForwardListError::InvalidKey) => (),
            Err(NestedForwardListError::KeyOutOfRange) => panic!(),
        }
    }
}
