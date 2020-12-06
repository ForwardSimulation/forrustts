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
    pub const fn new() -> NestedForwardList<Value> {
        return NestedForwardList {
            head_: Vec::<i32>::new(),
            tail_: Vec::<i32>::new(),
            next_: Vec::<i32>::new(),
            data_: Vec::<Value>::new(),
        };
    }

    pub fn null() -> i32 {
        return -1;
    }

    pub fn extend(&mut self, k: i32, v: Value) -> () {
        let idx = k as usize;
        if idx >= self.head_.len() {
            self.head_
                .resize(idx + 1, NestedForwardList::<Value>::null());
            self.tail_
                .resize(idx + 1, NestedForwardList::<Value>::null());
        }

        if self.head_[idx] == NestedForwardList::<Value>::null() {
            self.insert_new_record(idx as i32, v);
            return;
        }
        let t = self.tail_[idx];
        if t == NestedForwardList::<Value>::null() {
            // FIXME: need to error out here
        }
        self.data_.push(v);
        self.tail_[idx] = (self.data_.len() - 1) as i32;
        // NOTE: this differs from fwdpp, to avoid cast,
        // and thus could be failure point?
        self.next_[t as usize] = self.tail_[idx];
        self.next_.push(NestedForwardList::<Value>::null());
    }

    pub fn insert_new_record(&mut self, k: i32, v: Value) -> () {
        self.data_.push(v);
        let x = self.data_.len() as i32;
        self.head_[k as usize] = x;
        self.tail_[k as usize] = self.head_[k as usize];
        self.next_.push(NestedForwardList::<Value>::null());
    }

    pub fn fetch_mut(&mut self, at: i32) -> &mut Value {
        return &mut self.data_[at as usize];
    }

    pub fn fetch(&self, at: i32) -> &Value {
        return &self.data_[at as usize];
    }

    pub fn clear(&mut self) {
        self.data_.clear();
        self.head_.clear();
        self.tail_.clear();
        self.next_.clear();
    }

    pub fn nullify_list(&mut self, at: i32) -> () {
        self.head_[at as usize] = NestedForwardList::<Value>::null();
        self.tail_[at as usize] = NestedForwardList::<Value>::null();
    }

    pub fn reset(&mut self, newsize: usize) -> () {
        self.clear();
        self.head_
            .resize(newsize, NestedForwardList::<Value>::null());
        self.tail_
            .resize(newsize, NestedForwardList::<Value>::null());
        // .fill() is experimental right now...
        //self.head_.fill(NestedForwardList::<Value>::null());
        //self.tail_.fill(NestedForwardList::<Value>::null());
        // ... so we do this lambda mapping instead
        self.head_
            .iter_mut()
            .map(|x| *x = NestedForwardList::<Value>::null())
            .count();
        self.tail_
            .iter_mut()
            .map(|x| *x = NestedForwardList::<Value>::null())
            .count();
    }
}
