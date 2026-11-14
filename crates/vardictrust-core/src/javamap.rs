use std::fmt::Debug;
use std::hash::Hash;

const DEFAULT_CAPACITY: usize = 16;
const TREEIFY_THRESHOLD: usize = 8;
const MIN_TREEIFY_CAPACITY: usize = 64;
const DEFAULT_LOAD_FACTOR: f32 = 0.75;
const NONE: u32 = u32::MAX;

#[derive(Clone, Debug)]
struct Node<K, V> {
    key: K,
    val: V,
    next: u32,
}

#[derive(Debug)]
struct NewKeyLog<K>(Option<Vec<K>>);

impl<K> Default for NewKeyLog<K> {
    fn default() -> Self {
        Self(None)
    }
}

impl<K> Clone for NewKeyLog<K> {
    fn clone(&self) -> Self {
        Self(None)
    }
}

#[derive(Clone, Debug)]
pub struct JavaHashMap<K, V> {
    buckets: Vec<u32>,
    nodes: Vec<Node<K, V>>,
    free_head: u32,
    size: usize,
    threshold: usize,
    treeified: bool,
    fallback: V,
    new_key_log: NewKeyLog<K>,
}

#[derive(Debug)]
pub struct Iter<'a, K, V> {
    buckets: &'a [u32],
    nodes: &'a [Node<K, V>],
    bucket_idx: usize,
    curr: u32,
}

impl<'a, K, V> Iterator for Iter<'a, K, V> {
    type Item = (&'a K, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        while self.curr == NONE {
            let &head = self.buckets.get(self.bucket_idx)?;
            self.bucket_idx += 1;
            self.curr = head;
        }
        let node = self.nodes.get(self.curr as usize)?;
        self.curr = node.next;
        Some((&node.key, &node.val))
    }
}

/// Computes the 32-bit signed `Object.hashCode()` value used by [`JavaHashMap`].
pub trait JavaHash {
    fn java_hash(&self) -> i32;
}

impl JavaHash for i32 {
    #[inline]
    fn java_hash(&self) -> i32 {
        *self
    }
}

impl JavaHash for str {
    #[inline]
    fn java_hash(&self) -> i32 {
        java_string_hash(self)
    }
}

impl JavaHash for String {
    #[inline]
    fn java_hash(&self) -> i32 {
        java_string_hash(self)
    }
}

impl<T: JavaHash + ?Sized> JavaHash for &T {
    #[inline]
    fn java_hash(&self) -> i32 {
        (**self).java_hash()
    }
}

impl<K, V> JavaHashMap<K, V>
where
    K: Copy + JavaHash + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    pub fn new() -> Self {
        Self {
            buckets: Vec::new(),
            nodes: Vec::new(),
            free_head: NONE,
            size: 0,
            threshold: 0,
            treeified: false,
            fallback: V::default(),
            new_key_log: NewKeyLog::default(),
        }
    }

    /// Records newly allocated keys during an insert-only window (`sv::partial_pipeline`).
    /// Idempotent when already enabled; `remove` within the window is outside the contract.
    pub(crate) fn enable_new_key_log(&mut self) {
        self.new_key_log.0.get_or_insert_with(Vec::new);
    }

    pub(crate) fn take_new_key_log(&mut self) -> Vec<K> {
        self.new_key_log.0.take().unwrap_or_default()
    }

    fn spread(key: K) -> u32 {
        let h = key.java_hash() as u32;
        h ^ (h >> 16)
    }

    fn alloc_node(&mut self, key: K, val: V, next: u32) -> u32 {
        if let Some(ref mut log) = self.new_key_log.0 {
            log.push(key);
        }
        if self.free_head != NONE {
            let idx = self.free_head;
            if let Some(slot) = self.nodes.get_mut(idx as usize) {
                self.free_head = slot.next;
                *slot = Node { key, val, next };
                return idx;
            }
            self.free_head = NONE;
        }
        let idx = self.nodes.len() as u32;
        self.nodes.push(Node { key, val, next });
        idx
    }

    pub fn treeified(&self) -> bool {
        self.treeified
    }

    pub fn capacity(&self) -> usize {
        self.buckets.len()
    }

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// Pre-allocates capacity in the backing `nodes` arena without altering bucket
    /// count, resize threshold, or Java `HashMap` bucket-chain iteration order (`E4`).
    pub fn reserve_nodes(&mut self, additional: usize) {
        self.nodes.reserve(additional);
    }

    pub fn contains_key(&self, key: &K) -> bool {
        self.get(key).is_some()
    }

    /// Returns `true` iff both `first` and `second` are present in the map and `first`
    /// strictly precedes `second` in [`Self::iter`] order (`false` when `first == second`).
    pub fn contains_both_and_first_precedes(&self, first: &K, second: &K) -> bool {
        if first == second || self.buckets.is_empty() {
            return false;
        }
        let n = self.buckets.len();
        let idx_first = (Self::spread(*first) as usize) & (n - 1);
        let idx_second = (Self::spread(*second) as usize) & (n - 1);
        match idx_first.cmp(&idx_second) {
            std::cmp::Ordering::Greater => false,
            std::cmp::Ordering::Less => self.contains_key(first) && self.contains_key(second),
            std::cmp::Ordering::Equal => {
                let mut curr = self.buckets.get(idx_first).copied().unwrap_or(NONE);
                let mut saw_first = false;
                while curr != NONE {
                    let Some(node) = self.nodes.get(curr as usize) else {
                        return false;
                    };
                    if &node.key == first {
                        saw_first = true;
                    } else if &node.key == second {
                        return saw_first;
                    }
                    curr = node.next;
                }
                false
            }
        }
    }

    pub fn get(&self, key: &K) -> Option<&V> {
        if self.buckets.is_empty() {
            return None;
        }
        let hash = Self::spread(*key);
        let n = self.buckets.len();
        let idx = (hash as usize) & (n - 1);
        let mut curr = self.buckets.get(idx).copied().unwrap_or(NONE);
        while curr != NONE {
            let node = self.nodes.get(curr as usize)?;
            if &node.key == key {
                return Some(&node.val);
            }
            curr = node.next;
        }
        None
    }

    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        if self.buckets.is_empty() {
            return None;
        }
        let hash = Self::spread(*key);
        let n = self.buckets.len();
        let idx = (hash as usize) & (n - 1);
        let mut curr = self.buckets.get(idx).copied().unwrap_or(NONE);
        while curr != NONE {
            let node = self.nodes.get(curr as usize)?;
            if &node.key == key {
                return self.nodes.get_mut(curr as usize).map(|n| &mut n.val);
            }
            curr = node.next;
        }
        None
    }

    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if self.buckets.is_empty() {
            self.resize();
        }
        let hash = Self::spread(key);
        let n = self.buckets.len();
        let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };

        let mut curr = self.buckets.get(idx).copied().unwrap_or(NONE);
        let mut prev = NONE;
        let mut bin_len = 0usize;
        while curr != NONE {
            let Some(node) = self.nodes.get(curr as usize) else {
                break;
            };
            if node.key == key {
                if let Some(node_mut) = self.nodes.get_mut(curr as usize) {
                    let old = std::mem::replace(&mut node_mut.val, value);
                    return Some(old);
                }
                return None;
            }
            prev = curr;
            curr = node.next;
            bin_len += 1;
        }

        let new_idx = self.alloc_node(key, value, NONE);
        if prev == NONE {
            if let Some(b) = self.buckets.get_mut(idx) {
                *b = new_idx;
            }
        } else if let Some(prev_node) = self.nodes.get_mut(prev as usize) {
            prev_node.next = new_idx;
        }

        if bin_len >= TREEIFY_THRESHOLD {
            self.treeify_bin();
        }
        self.size += 1;
        if self.size > self.threshold {
            self.resize();
        }
        None
    }

    /// Emulates `V v = m.get(k); if (v == null) { v = f(); m.put(k, v); }` (`HashMap.putVal`):
    /// single bucket walk, tail append, treeify at >= 8 existing nodes (9th node),
    /// and resize after insert when `size > threshold`.
    pub fn get_or_put_with<F>(&mut self, key: K, f: F) -> &mut V
    where
        F: FnOnce() -> V,
    {
        if self.buckets.is_empty() {
            self.resize();
        }
        let hash = Self::spread(key);
        let n = self.buckets.len();
        let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };

        let mut curr = self.buckets.get(idx).copied().unwrap_or(NONE);
        let mut prev = NONE;
        let mut bin_len = 0usize;
        while curr != NONE {
            let Some(node) = self.nodes.get(curr as usize) else {
                break;
            };
            if node.key == key {
                return self
                    .nodes
                    .get_mut(curr as usize)
                    .map(|n| &mut n.val)
                    .unwrap_or(&mut self.fallback);
            }
            prev = curr;
            curr = node.next;
            bin_len += 1;
        }

        let new_idx = self.alloc_node(key, f(), NONE);
        if prev == NONE {
            if let Some(b) = self.buckets.get_mut(idx) {
                *b = new_idx;
            }
        } else if let Some(prev_node) = self.nodes.get_mut(prev as usize) {
            prev_node.next = new_idx;
        }

        if bin_len >= TREEIFY_THRESHOLD {
            self.treeify_bin();
        }
        self.size += 1;
        if self.size > self.threshold {
            self.resize();
        }

        if let Some(n) = self.nodes.get_mut(new_idx as usize) {
            &mut n.val
        } else {
            &mut self.fallback
        }
    }

    pub fn compute_if_absent<F>(&mut self, key: K, f: F) -> &mut V
    where
        F: FnOnce(&K) -> V,
    {
        if self.buckets.is_empty() || self.size > self.threshold {
            // OpenJDK HashMap.computeIfAbsent pre-resize check (size > threshold)
            self.resize();
        }
        let hash = Self::spread(key);
        let n = self.buckets.len();
        let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };

        let mut curr = self.buckets.get(idx).copied().unwrap_or(NONE);
        let mut bin_len = 0usize;
        let mut found = NONE;
        while curr != NONE {
            let Some(node) = self.nodes.get(curr as usize) else {
                break;
            };
            if node.key == key {
                found = curr;
                break;
            }
            curr = node.next;
            bin_len += 1;
        }

        let target_idx = if found != NONE {
            found
        } else {
            let val = f(&key);
            let head = self.buckets.get(idx).copied().unwrap_or(NONE);
            // OpenJDK HashMap.computeIfAbsent inserts at head of bin
            let new_idx = self.alloc_node(key, val, head);
            if let Some(b) = self.buckets.get_mut(idx) {
                *b = new_idx;
            }
            self.size += 1;
            if bin_len + 1 >= TREEIFY_THRESHOLD {
                self.treeify_bin();
            }
            new_idx
        };

        if let Some(n) = self.nodes.get_mut(target_idx as usize) {
            &mut n.val
        } else {
            &mut self.fallback
        }
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
        if self.buckets.is_empty() {
            return None;
        }
        let hash = Self::spread(*key);
        let n = self.buckets.len();
        let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };

        let mut curr = self.buckets.get(idx).copied().unwrap_or(NONE);
        let mut prev = NONE;
        while curr != NONE {
            let Some(node) = self.nodes.get(curr as usize) else {
                break;
            };
            if &node.key == key {
                let next = node.next;
                if prev == NONE {
                    if let Some(b) = self.buckets.get_mut(idx) {
                        *b = next;
                    }
                } else if let Some(prev_node) = self.nodes.get_mut(prev as usize) {
                    prev_node.next = next;
                }
                self.size = self.size.saturating_sub(1);
                if let Some(removed_node) = self.nodes.get_mut(curr as usize) {
                    removed_node.next = self.free_head;
                    self.free_head = curr;
                    return Some(std::mem::take(&mut removed_node.val));
                }
                return None;
            }
            prev = curr;
            curr = node.next;
        }
        None
    }

    fn treeify_bin(&mut self) {
        if self.buckets.len() < MIN_TREEIFY_CAPACITY {
            self.resize();
        } else {
            self.treeified = true;
        }
    }

    fn resize(&mut self) {
        let old_cap = self.buckets.len();
        if old_cap == 0 {
            self.buckets = vec![NONE; DEFAULT_CAPACITY];
            self.threshold = (DEFAULT_CAPACITY as f32 * DEFAULT_LOAD_FACTOR) as usize;
            self.nodes.reserve(DEFAULT_CAPACITY);
            return;
        }

        let new_cap = old_cap << 1;
        self.threshold = (new_cap as f32 * DEFAULT_LOAD_FACTOR) as usize;
        self.buckets.resize(new_cap, NONE);

        for j in 0..old_cap {
            let mut curr = self.buckets.get(j).copied().unwrap_or(NONE);
            if curr == NONE {
                continue;
            }
            let mut lo_head = NONE;
            let mut lo_tail = NONE;
            let mut hi_head = NONE;
            let mut hi_tail = NONE;

            while curr != NONE {
                let Some(node) = self.nodes.get(curr as usize) else {
                    break;
                };
                let next = node.next;
                let hash = Self::spread(node.key);
                if (hash & (old_cap as u32)) == 0 {
                    if lo_tail == NONE {
                        lo_head = curr;
                    } else if let Some(prev_node) = self.nodes.get_mut(lo_tail as usize) {
                        prev_node.next = curr;
                    }
                    lo_tail = curr;
                } else {
                    if hi_tail == NONE {
                        hi_head = curr;
                    } else if let Some(prev_node) = self.nodes.get_mut(hi_tail as usize) {
                        prev_node.next = curr;
                    }
                    hi_tail = curr;
                }
                curr = next;
            }

            if lo_tail != NONE {
                if let Some(tail_node) = self.nodes.get_mut(lo_tail as usize) {
                    tail_node.next = NONE;
                }
            }
            if let Some(b) = self.buckets.get_mut(j) {
                *b = lo_head;
            }

            if hi_tail != NONE {
                if let Some(tail_node) = self.nodes.get_mut(hi_tail as usize) {
                    tail_node.next = NONE;
                }
            }
            if let Some(b) = self.buckets.get_mut(j + old_cap) {
                *b = hi_head;
            }
        }
    }

    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.iter().map(|(k, _)| k)
    }

    pub fn iter(&self) -> Iter<'_, K, V> {
        Iter {
            buckets: &self.buckets,
            nodes: &self.nodes,
            bucket_idx: 0,
            curr: NONE,
        }
    }
}

impl<K, V> Default for JavaHashMap<K, V>
where
    K: Copy + JavaHash + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<K, V> IntoIterator for JavaHashMap<K, V>
where
    K: Copy + JavaHash + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    type Item = (K, V);
    type IntoIter = std::vec::IntoIter<(K, V)>;

    fn into_iter(mut self) -> Self::IntoIter {
        let mut items = Vec::with_capacity(self.size);
        for &head in &self.buckets {
            let mut curr = head;
            while curr != NONE {
                let Some(node) = self.nodes.get_mut(curr as usize) else {
                    break;
                };
                let k = node.key;
                let v = std::mem::take(&mut node.val);
                curr = node.next;
                items.push((k, v));
            }
        }
        items.into_iter()
    }
}

/// Computes `java.lang.String.hashCode()` over UTF-16 code units:
/// `s[0]*31^(n-1) + s[1]*31^(n-2) + ... + s[n-1]` using wrapping 32-bit signed arithmetic.
pub fn java_string_hash(s: &str) -> i32 {
    s.encode_utf16()
        .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(i32::from(c)))
}

/// Orders unique keys (given in first-insertion order) in OpenJDK 8+ `HashMap<String, V>`
/// `keySet()` / `entrySet()` iteration order, matching `HashMap.putVal` bucket chaining,
/// `treeifyBin` capacity-doubling when a bucket reaches 9 entries (`cap < 64`),
/// `TreeNode.treeify` / `putTreeVal` / `moveRootToFront` red-black tree ordering (`cap >= 64`),
/// and `loHead`/`hiHead` (`TreeNode.split`) on resize.
pub fn java_hashmap_key_order<'a, S, I>(keys_in_insertion_order: I) -> Vec<&'a S>
where
    S: AsRef<str> + 'a + ?Sized,
    I: IntoIterator<Item = &'a S>,
{
    const UNTREEIFY_THRESHOLD: usize = 6;

    struct StrNode<'a, S: ?Sized> {
        key: &'a S,
        hash: i32,
        next: u32,
        prev: u32,
        parent: u32,
        left: u32,
        right: u32,
        red: bool,
        is_tree: bool,
    }

    fn spread_str(s: &str) -> i32 {
        let h = java_string_hash(s);
        h ^ ((h as u32 >> 16) as i32)
    }

    fn cmp_str(h_k: i32, k: &str, h_p: i32, p: &str) -> i32 {
        if h_p > h_k {
            -1
        } else if h_p < h_k {
            1
        } else {
            match k.encode_utf16().cmp(p.encode_utf16()) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            }
        }
    }

    fn move_root_to_front<S: ?Sized>(
        buckets: &mut [u32],
        nodes: &mut [StrNode<'_, S>],
        bucket_idx: usize,
        root: u32,
    ) {
        let first = buckets.get(bucket_idx).copied().unwrap_or(NONE);
        if root != NONE && root != first {
            if let Some(b) = buckets.get_mut(bucket_idx) {
                *b = root;
            }
            let (rp, rn) = nodes
                .get(root as usize)
                .map_or((NONE, NONE), |n| (n.prev, n.next));
            if rn != NONE {
                if let Some(rn_node) = nodes.get_mut(rn as usize) {
                    rn_node.prev = rp;
                }
            }
            if rp != NONE {
                if let Some(rp_node) = nodes.get_mut(rp as usize) {
                    rp_node.next = rn;
                }
            }
            if first != NONE {
                if let Some(first_node) = nodes.get_mut(first as usize) {
                    first_node.prev = root;
                }
            }
            if let Some(root_node) = nodes.get_mut(root as usize) {
                root_node.next = first;
                root_node.prev = NONE;
            }
        }
    }

    fn rotate_left<S: ?Sized>(nodes: &mut [StrNode<'_, S>], mut root: u32, p: u32) -> u32 {
        if p != NONE {
            let r = nodes.get(p as usize).map_or(NONE, |n| n.right);
            if r != NONE {
                let rl = nodes.get(r as usize).map_or(NONE, |n| n.left);
                if let Some(pn) = nodes.get_mut(p as usize) {
                    pn.right = rl;
                }
                if rl != NONE {
                    if let Some(rln) = nodes.get_mut(rl as usize) {
                        rln.parent = p;
                    }
                }
                let pp = nodes.get(p as usize).map_or(NONE, |n| n.parent);
                if let Some(rn) = nodes.get_mut(r as usize) {
                    rn.parent = pp;
                }
                if pp == NONE {
                    root = r;
                    if let Some(rn) = nodes.get_mut(r as usize) {
                        rn.red = false;
                    }
                } else if nodes.get(pp as usize).is_some_and(|n| n.left == p) {
                    if let Some(ppn) = nodes.get_mut(pp as usize) {
                        ppn.left = r;
                    }
                } else if let Some(ppn) = nodes.get_mut(pp as usize) {
                    ppn.right = r;
                }
                if let Some(rn) = nodes.get_mut(r as usize) {
                    rn.left = p;
                }
                if let Some(pn) = nodes.get_mut(p as usize) {
                    pn.parent = r;
                }
            }
        }
        root
    }

    fn rotate_right<S: ?Sized>(nodes: &mut [StrNode<'_, S>], mut root: u32, p: u32) -> u32 {
        if p != NONE {
            let l = nodes.get(p as usize).map_or(NONE, |n| n.left);
            if l != NONE {
                let lr = nodes.get(l as usize).map_or(NONE, |n| n.right);
                if let Some(pn) = nodes.get_mut(p as usize) {
                    pn.left = lr;
                }
                if lr != NONE {
                    if let Some(lrn) = nodes.get_mut(lr as usize) {
                        lrn.parent = p;
                    }
                }
                let pp = nodes.get(p as usize).map_or(NONE, |n| n.parent);
                if let Some(ln) = nodes.get_mut(l as usize) {
                    ln.parent = pp;
                }
                if pp == NONE {
                    root = l;
                    if let Some(ln) = nodes.get_mut(l as usize) {
                        ln.red = false;
                    }
                } else if nodes.get(pp as usize).is_some_and(|n| n.right == p) {
                    if let Some(ppn) = nodes.get_mut(pp as usize) {
                        ppn.right = l;
                    }
                } else if let Some(ppn) = nodes.get_mut(pp as usize) {
                    ppn.left = l;
                }
                if let Some(ln) = nodes.get_mut(l as usize) {
                    ln.right = p;
                }
                if let Some(pn) = nodes.get_mut(p as usize) {
                    pn.parent = l;
                }
            }
        }
        root
    }

    fn balance_insertion<S: ?Sized>(
        nodes: &mut [StrNode<'_, S>],
        mut root: u32,
        mut x: u32,
    ) -> u32 {
        if let Some(xn) = nodes.get_mut(x as usize) {
            xn.red = true;
        }
        loop {
            let xp = nodes.get(x as usize).map_or(NONE, |n| n.parent);
            if xp == NONE {
                if let Some(xn) = nodes.get_mut(x as usize) {
                    xn.red = false;
                }
                return x;
            }
            let (xp_red, xpp) = nodes
                .get(xp as usize)
                .map_or((false, NONE), |n| (n.red, n.parent));
            if !xp_red || xpp == NONE {
                return root;
            }
            let (xppl, xppr) = nodes
                .get(xpp as usize)
                .map_or((NONE, NONE), |n| (n.left, n.right));
            if xp == xppl {
                if xppr != NONE && nodes.get(xppr as usize).is_some_and(|n| n.red) {
                    if let Some(n) = nodes.get_mut(xppr as usize) {
                        n.red = false;
                    }
                    if let Some(n) = nodes.get_mut(xp as usize) {
                        n.red = false;
                    }
                    if let Some(n) = nodes.get_mut(xpp as usize) {
                        n.red = true;
                    }
                    x = xpp;
                } else {
                    if nodes.get(xp as usize).is_some_and(|n| x == n.right) {
                        x = xp;
                        root = rotate_left(nodes, root, x);
                    }
                    let xp = nodes.get(x as usize).map_or(NONE, |n| n.parent);
                    let xpp = nodes.get(xp as usize).map_or(NONE, |n| n.parent);
                    if xp != NONE {
                        if let Some(n) = nodes.get_mut(xp as usize) {
                            n.red = false;
                        }
                        if xpp != NONE {
                            if let Some(n) = nodes.get_mut(xpp as usize) {
                                n.red = true;
                            }
                            root = rotate_right(nodes, root, xpp);
                        }
                    }
                }
            } else if xppl != NONE && nodes.get(xppl as usize).is_some_and(|n| n.red) {
                if let Some(n) = nodes.get_mut(xppl as usize) {
                    n.red = false;
                }
                if let Some(n) = nodes.get_mut(xp as usize) {
                    n.red = false;
                }
                if let Some(n) = nodes.get_mut(xpp as usize) {
                    n.red = true;
                }
                x = xpp;
            } else {
                if nodes.get(xp as usize).is_some_and(|n| x == n.left) {
                    x = xp;
                    root = rotate_right(nodes, root, x);
                }
                let xp = nodes.get(x as usize).map_or(NONE, |n| n.parent);
                let xpp = nodes.get(xp as usize).map_or(NONE, |n| n.parent);
                if xp != NONE {
                    if let Some(n) = nodes.get_mut(xp as usize) {
                        n.red = false;
                    }
                    if xpp != NONE {
                        if let Some(n) = nodes.get_mut(xpp as usize) {
                            n.red = true;
                        }
                        root = rotate_left(nodes, root, xpp);
                    }
                }
            }
        }
    }

    fn treeify_bucket<S: AsRef<str> + ?Sized>(
        buckets: &mut [u32],
        nodes: &mut [StrNode<'_, S>],
        bucket_idx: usize,
    ) {
        let head = buckets.get(bucket_idx).copied().unwrap_or(NONE);
        if head == NONE {
            return;
        }
        let mut prev = NONE;
        let mut curr = head;
        while curr != NONE {
            let next = if let Some(n) = nodes.get_mut(curr as usize) {
                n.prev = prev;
                n.parent = NONE;
                n.left = NONE;
                n.right = NONE;
                n.red = false;
                n.is_tree = true;
                n.next
            } else {
                break;
            };
            prev = curr;
            curr = next;
        }

        let mut root = NONE;
        let mut x = head;
        while x != NONE {
            let next = nodes.get(x as usize).map_or(NONE, |n| n.next);
            if root == NONE {
                if let Some(xn) = nodes.get_mut(x as usize) {
                    xn.parent = NONE;
                    xn.red = false;
                }
                root = x;
            } else {
                let (h, k_ref) = match nodes.get(x as usize) {
                    Some(n) => (n.hash, n.key),
                    None => break,
                };
                let mut p = root;
                while let Some(n) = nodes.get(p as usize) {
                    let (ph, pk_ref, pl, pr) = (n.hash, n.key, n.left, n.right);
                    let dir = cmp_str(h, k_ref.as_ref(), ph, pk_ref.as_ref());
                    let xp = p;
                    p = if dir <= 0 { pl } else { pr };
                    if p == NONE {
                        if let Some(xn) = nodes.get_mut(x as usize) {
                            xn.parent = xp;
                        }
                        if let Some(xpn) = nodes.get_mut(xp as usize) {
                            if dir <= 0 {
                                xpn.left = x;
                            } else {
                                xpn.right = x;
                            }
                        }
                        root = balance_insertion(nodes, root, x);
                        break;
                    }
                }
            }
            x = next;
        }
        move_root_to_front(buckets, nodes, bucket_idx, root);
    }

    fn resize<S: AsRef<str> + ?Sized>(
        buckets: &mut Vec<u32>,
        nodes: &mut [StrNode<'_, S>],
        threshold: &mut usize,
    ) {
        let old_cap = buckets.len();
        if old_cap == 0 {
            *buckets = vec![NONE; DEFAULT_CAPACITY];
            *threshold = (DEFAULT_CAPACITY as f32 * DEFAULT_LOAD_FACTOR) as usize;
            return;
        }
        let new_cap = old_cap << 1;
        *threshold = (new_cap as f32 * DEFAULT_LOAD_FACTOR) as usize;
        buckets.resize(new_cap, NONE);

        for j in 0..old_cap {
            let e = buckets.get(j).copied().unwrap_or(NONE);
            if e == NONE {
                continue;
            }
            if let Some(b) = buckets.get_mut(j) {
                *b = NONE;
            }
            let is_tree = nodes.get(e as usize).is_some_and(|n| n.is_tree);
            let mut lo_head = NONE;
            let mut lo_tail = NONE;
            let mut hi_head = NONE;
            let mut hi_tail = NONE;
            let mut lc = 0usize;
            let mut hc = 0usize;
            let mut curr = e;
            while curr != NONE {
                let (next, hash) = match nodes.get(curr as usize) {
                    Some(n) => (n.next, n.hash as u32),
                    None => break,
                };
                if (hash & (old_cap as u32)) == 0 {
                    if let Some(cn) = nodes.get_mut(curr as usize) {
                        cn.prev = lo_tail;
                    }
                    if lo_tail == NONE {
                        lo_head = curr;
                    } else if let Some(lt) = nodes.get_mut(lo_tail as usize) {
                        lt.next = curr;
                    }
                    lo_tail = curr;
                    lc += 1;
                } else {
                    if let Some(cn) = nodes.get_mut(curr as usize) {
                        cn.prev = hi_tail;
                    }
                    if hi_tail == NONE {
                        hi_head = curr;
                    } else if let Some(ht) = nodes.get_mut(hi_tail as usize) {
                        ht.next = curr;
                    }
                    hi_tail = curr;
                    hc += 1;
                }
                curr = next;
            }
            if lo_tail != NONE {
                if let Some(lt) = nodes.get_mut(lo_tail as usize) {
                    lt.next = NONE;
                }
                if let Some(b) = buckets.get_mut(j) {
                    *b = lo_head;
                }
                if is_tree {
                    if lc <= UNTREEIFY_THRESHOLD {
                        let mut c = lo_head;
                        while c != NONE {
                            let next = if let Some(n) = nodes.get_mut(c as usize) {
                                n.is_tree = false;
                                n.next
                            } else {
                                break;
                            };
                            c = next;
                        }
                    } else if hi_head != NONE {
                        treeify_bucket(buckets, nodes, j);
                    }
                }
            }
            if hi_tail != NONE {
                if let Some(ht) = nodes.get_mut(hi_tail as usize) {
                    ht.next = NONE;
                }
                if let Some(b) = buckets.get_mut(j + old_cap) {
                    *b = hi_head;
                }
                if is_tree {
                    if hc <= UNTREEIFY_THRESHOLD {
                        let mut c = hi_head;
                        while c != NONE {
                            let next = if let Some(n) = nodes.get_mut(c as usize) {
                                n.is_tree = false;
                                n.next
                            } else {
                                break;
                            };
                            c = next;
                        }
                    } else if lo_head != NONE {
                        treeify_bucket(buckets, nodes, j + old_cap);
                    }
                }
            }
        }
    }

    let mut buckets: Vec<u32> = Vec::new();
    let mut nodes: Vec<StrNode<'a, S>> = Vec::new();
    let mut size = 0usize;
    let mut threshold = 0usize;

    for key in keys_in_insertion_order {
        if buckets.is_empty() {
            resize(&mut buckets, &mut nodes, &mut threshold);
        }
        let hash = spread_str(key.as_ref());
        let n = buckets.len();
        let idx = (hash as u32 as usize) & (n - 1);
        let first = buckets.get(idx).copied().unwrap_or(NONE);

        if first == NONE {
            let new_idx = nodes.len() as u32;
            nodes.push(StrNode {
                key,
                hash,
                next: NONE,
                prev: NONE,
                parent: NONE,
                left: NONE,
                right: NONE,
                red: false,
                is_tree: false,
            });
            if let Some(b) = buckets.get_mut(idx) {
                *b = new_idx;
            }
        } else if nodes.get(first as usize).is_some_and(|nd| nd.is_tree) {
            let mut p = first;
            let mut dup = false;
            let mut insert_at = (NONE, 0i32);
            while let Some(nd) = nodes.get(p as usize) {
                let (ph, pk_ref, pl, pr) = (nd.hash, nd.key, nd.left, nd.right);
                if ph == hash && pk_ref.as_ref() == key.as_ref() {
                    dup = true;
                    break;
                }
                let dir = cmp_str(hash, key.as_ref(), ph, pk_ref.as_ref());
                let xp = p;
                p = if dir <= 0 { pl } else { pr };
                if p == NONE {
                    insert_at = (xp, dir);
                    break;
                }
            }
            if dup {
                continue;
            }
            let (xp, dir) = insert_at;
            let xpn = nodes.get(xp as usize).map_or(NONE, |nd| nd.next);
            let new_idx = nodes.len() as u32;
            nodes.push(StrNode {
                key,
                hash,
                next: xpn,
                prev: xp,
                parent: xp,
                left: NONE,
                right: NONE,
                red: true,
                is_tree: true,
            });
            if let Some(xpn_node) = nodes.get_mut(xp as usize) {
                xpn_node.next = new_idx;
                if dir <= 0 {
                    xpn_node.left = new_idx;
                } else {
                    xpn_node.right = new_idx;
                }
            }
            if xpn != NONE {
                if let Some(next_node) = nodes.get_mut(xpn as usize) {
                    next_node.prev = new_idx;
                }
            }
            let root = balance_insertion(&mut nodes, first, new_idx);
            move_root_to_front(&mut buckets, &mut nodes, idx, root);
        } else {
            let mut curr = first;
            let mut prev = NONE;
            let mut bin_len = 0usize;
            let mut dup = false;
            while curr != NONE {
                let Some(nd) = nodes.get(curr as usize) else {
                    break;
                };
                if nd.hash == hash && nd.key.as_ref() == key.as_ref() {
                    dup = true;
                    break;
                }
                prev = curr;
                curr = nd.next;
                bin_len += 1;
            }
            if dup {
                continue;
            }
            let new_idx = nodes.len() as u32;
            nodes.push(StrNode {
                key,
                hash,
                next: NONE,
                prev,
                parent: NONE,
                left: NONE,
                right: NONE,
                red: false,
                is_tree: false,
            });
            if let Some(pn) = nodes.get_mut(prev as usize) {
                pn.next = new_idx;
            }
            if bin_len >= TREEIFY_THRESHOLD {
                if buckets.len() < MIN_TREEIFY_CAPACITY {
                    resize(&mut buckets, &mut nodes, &mut threshold);
                } else {
                    treeify_bucket(&mut buckets, &mut nodes, idx);
                }
            }
        }

        size += 1;
        if size > threshold {
            resize(&mut buckets, &mut nodes, &mut threshold);
        }
    }

    let mut out = Vec::with_capacity(size);
    for &head in &buckets {
        let mut curr = head;
        while curr != NONE {
            let Some(nd) = nodes.get(curr as usize) else {
                break;
            };
            out.push(nd.key);
            curr = nd.next;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Debug)]
    struct JavaHashMapLegacy<K, V> {
        table: Vec<Vec<(K, V)>>,
        size: usize,
        threshold: usize,
        treeified: bool,
    }

    impl<K, V> JavaHashMapLegacy<K, V>
    where
        K: Copy + Into<i32> + Eq + PartialEq + Hash + Debug,
        V: Clone + Debug + Default,
    {
        fn new() -> Self {
            Self {
                table: Vec::new(),
                size: 0,
                threshold: 0,
                treeified: false,
            }
        }

        fn spread(key: K) -> u32 {
            let h = key.into() as u32;
            h ^ (h >> 16)
        }

        fn treeified(&self) -> bool {
            self.treeified
        }

        fn len(&self) -> usize {
            self.size
        }

        fn is_empty(&self) -> bool {
            self.size == 0
        }

        fn contains_key(&self, key: &K) -> bool {
            self.get(key).is_some()
        }

        fn get(&self, key: &K) -> Option<&V> {
            if self.table.is_empty() {
                return None;
            }
            let hash = Self::spread(*key);
            let n = self.table.len();
            let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };
            if let Some(bin) = self.table.get(idx) {
                for (k, v) in bin {
                    if k == key {
                        return Some(v);
                    }
                }
            }
            None
        }

        fn get_mut(&mut self, key: &K) -> Option<&mut V> {
            if self.table.is_empty() {
                return None;
            }
            let hash = Self::spread(*key);
            let n = self.table.len();
            let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };
            if let Some(bin) = self.table.get_mut(idx) {
                for (k, v) in bin {
                    if k == key {
                        return Some(v);
                    }
                }
            }
            None
        }

        fn insert(&mut self, key: K, value: V) -> Option<V> {
            if self.table.is_empty() {
                self.resize();
            }
            let hash = Self::spread(key);
            let n = self.table.len();
            let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };

            if let Some(bin) = self.table.get_mut(idx) {
                for entry in bin.iter_mut() {
                    if entry.0 == key {
                        let old = std::mem::replace(&mut entry.1, value);
                        return Some(old);
                    }
                }

                bin.push((key, value));
                let bin_len = bin.len();
                if bin_len > TREEIFY_THRESHOLD {
                    self.treeify_bin();
                }
                self.size += 1;
                if self.size > self.threshold {
                    self.resize();
                }
            }
            None
        }

        fn get_or_put_with<F>(&mut self, key: K, f: F) -> &mut V
        where
            F: FnOnce() -> V,
        {
            if !self.contains_key(&key) {
                let val = f();
                self.insert(key, val);
            }
            self.get_mut(&key).expect("entry exists after insert")
        }

        fn compute_if_absent<F>(&mut self, key: K, f: F) -> &mut V
        where
            F: FnOnce(&K) -> V,
        {
            if self.table.is_empty() || self.size > self.threshold {
                self.resize();
            }
            let hash = Self::spread(key);
            let n = self.table.len();
            let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };

            let exists = match self.table.get(idx) {
                Some(bin) => bin.iter().any(|(k, _)| *k == key),
                None => false,
            };

            if !exists {
                let val = f(&key);
                if let Some(bin) = self.table.get_mut(idx) {
                    bin.insert(0, (key, val));
                    self.size += 1;
                    let bin_len = bin.len();
                    if bin_len >= TREEIFY_THRESHOLD {
                        self.treeify_bin();
                    }
                }
            }

            let hash = Self::spread(key);
            let n = self.table.len();
            let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };

            let pos = self
                .table
                .get(idx)
                .and_then(|bin| bin.iter().position(|(k, _)| *k == key));

            let (target_bin, target_pos) = match pos {
                Some(p) => (idx, p),
                None => {
                    let mut found = None;
                    for (b_idx, bin) in self.table.iter().enumerate() {
                        if let Some(p) = bin.iter().position(|(k, _)| *k == key) {
                            found = Some((b_idx, p));
                            break;
                        }
                    }
                    match found {
                        Some(fp) => fp,
                        None => (idx, 0),
                    }
                }
            };

            let valid = target_bin < self.table.len()
                && self
                    .table
                    .get(target_bin)
                    .is_some_and(|b| target_pos < b.len());
            if valid {
                let bin = &mut self.table[target_bin];
                return &mut bin[target_pos].1;
            }

            if self.table.is_empty() {
                self.table.push(Vec::new());
            }
            let bin = &mut self.table[0];
            if bin.is_empty() {
                bin.push((key, V::default()));
            }
            &mut bin[0].1
        }

        fn remove(&mut self, key: &K) -> Option<V> {
            if self.table.is_empty() {
                return None;
            }
            let hash = Self::spread(*key);
            let n = self.table.len();
            let idx = if n == 0 { 0 } else { (hash as usize) & (n - 1) };
            if let Some(bin) = self.table.get_mut(idx) {
                if let Some(pos) = bin.iter().position(|(k, _)| k == key) {
                    self.size = self.size.saturating_sub(1);
                    let (_, v) = bin.remove(pos);
                    return Some(v);
                }
            }
            None
        }

        fn treeify_bin(&mut self) {
            if self.table.len() < MIN_TREEIFY_CAPACITY {
                self.resize();
            } else {
                self.treeified = true;
            }
        }

        fn resize(&mut self) {
            let old_cap = self.table.len();
            let new_cap = if old_cap == 0 {
                DEFAULT_CAPACITY
            } else {
                old_cap << 1
            };

            self.threshold = (new_cap as f32 * DEFAULT_LOAD_FACTOR) as usize;

            if old_cap == 0 {
                self.table = vec![Vec::new(); new_cap];
                return;
            }

            let mut new_table = vec![Vec::new(); new_cap];
            for j in 0..old_cap {
                let mut lo = Vec::new();
                let mut hi = Vec::new();
                if let Some(bin) = self.table.get_mut(j) {
                    for (k, v) in std::mem::take(bin) {
                        let hash = Self::spread(k);
                        if (hash & (old_cap as u32)) == 0 {
                            lo.push((k, v));
                        } else {
                            hi.push((k, v));
                        }
                    }
                }
                if let Some(b) = new_table.get_mut(j) {
                    *b = lo;
                }
                if let Some(b) = new_table.get_mut(j + old_cap) {
                    *b = hi;
                }
            }
            self.table = new_table;
        }

        fn keys(&self) -> impl Iterator<Item = &K> {
            self.table.iter().flat_map(|bin| bin.iter().map(|(k, _)| k))
        }

        fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
            self.table
                .iter()
                .flat_map(|bin| bin.iter().map(|(k, v)| (k, v)))
        }

        fn into_vec(self) -> Vec<(K, V)> {
            let mut items = Vec::new();
            for mut bin in self.table {
                items.append(&mut bin);
            }
            items
        }
    }

    #[test]
    fn javamap_arena_matches_legacy_on_random_op_streams() {
        fn next_u64(state: &mut u64) -> u64 {
            let mut x = *state;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            *state = x;
            x
        }

        for stream_idx in 0..50u64 {
            let mut rng =
                0x9E37_79B9_7F4A_7C15u64 ^ (stream_idx.wrapping_mul(0xBF58_476D_1CE4_E5B9));
            if rng == 0 {
                rng = 1;
            }
            let mut arena: JavaHashMap<i32, i32> = JavaHashMap::new();
            let mut legacy: JavaHashMapLegacy<i32, i32> = JavaHashMapLegacy::new();

            for step in 0..5_000usize {
                let r = next_u64(&mut rng);
                let op = r % 100;
                // Mix general keys and heavy collision keys sharing low 16 bits
                // Note: spread(k) = k ^ (k >> 16). To force spread(k) to collide in low bits,
                // if k = (i << 16) | (i ^ 0x1234), then k ^ (k >> 16) has low 16 bits == 0x1234!
                // Also include raw `(i << 16) | 0x1234` keys as specified.
                let i = ((r >> 8) & 0x7F) as i32;
                let key = match (r >> 16) % 3 {
                    0 => (i << 16) | 0x1234,
                    1 => (i << 16) | ((i ^ 0x1234) & 0xFFFF),
                    _ => ((r >> 20) % 256) as i32 - 64,
                };
                let val = (step as i32).wrapping_mul(3) ^ key;

                if op < 30 {
                    let a_old = arena.insert(key, val);
                    let l_old = legacy.insert(key, val);
                    assert_eq!(
                        a_old, l_old,
                        "insert return mismatch at stream {stream_idx} step {step}"
                    );
                } else if op < 60 {
                    let delta = ((r >> 32) & 0xFF) as i32;
                    let a_ref = arena.compute_if_absent(key, |k| k.wrapping_add(val));
                    *a_ref = a_ref.wrapping_add(delta);
                    let l_ref = legacy.compute_if_absent(key, |k| k.wrapping_add(val));
                    *l_ref = l_ref.wrapping_add(delta);
                    assert_eq!(
                        *a_ref, *l_ref,
                        "compute_if_absent value mismatch at stream {stream_idx} step {step}"
                    );
                } else if op < 80 {
                    let a_rem = arena.remove(&key);
                    let l_rem = legacy.remove(&key);
                    assert_eq!(
                        a_rem, l_rem,
                        "remove mismatch at stream {stream_idx} step {step}"
                    );
                } else if op < 90 {
                    assert_eq!(arena.contains_key(&key), legacy.contains_key(&key));
                    assert_eq!(arena.get(&key), legacy.get(&key));
                } else if op < 97 {
                    let a_mut = arena.get_mut(&key);
                    let l_mut = legacy.get_mut(&key);
                    assert_eq!(a_mut.is_some(), l_mut.is_some());
                    if let (Some(av), Some(lv)) = (a_mut, l_mut) {
                        *av = av.wrapping_add(7);
                        *lv = lv.wrapping_add(7);
                        assert_eq!(*av, *lv);
                    }
                } else {
                    arena = arena.clone();
                    legacy = legacy.clone();
                }

                if step % 500 == 0 || step == 4_999 {
                    assert_eq!(
                        arena.len(),
                        legacy.len(),
                        "len mismatch at stream {stream_idx} step {step}"
                    );
                    assert_eq!(
                        arena.is_empty(),
                        legacy.is_empty(),
                        "is_empty mismatch at stream {stream_idx} step {step}"
                    );
                    assert_eq!(
                        arena.treeified(),
                        legacy.treeified(),
                        "treeified mismatch at stream {stream_idx} step {step}"
                    );
                    let a_keys: Vec<i32> = arena.keys().copied().collect();
                    let l_keys: Vec<i32> = legacy.keys().copied().collect();
                    assert_eq!(
                        a_keys, l_keys,
                        "keys() order mismatch at stream {stream_idx} step {step}"
                    );
                    let a_iter: Vec<(i32, i32)> = arena.iter().map(|(&k, &v)| (k, v)).collect();
                    let l_iter: Vec<(i32, i32)> = legacy.iter().map(|(&k, &v)| (k, v)).collect();
                    assert_eq!(
                        a_iter, l_iter,
                        "iter() mismatch at stream {stream_idx} step {step}"
                    );
                    let a_into: Vec<(i32, i32)> = arena.clone().into_iter().collect();
                    let l_into: Vec<(i32, i32)> = legacy.clone().into_vec();
                    assert_eq!(
                        a_into, l_into,
                        "into_iter() mismatch at stream {stream_idx} step {step}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_javamap_oracle() {
        let run_oracle = |script: &str, expected_keys: Option<&str>| {
            let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
            for line in script.lines() {
                if let Some(k_str) = line.strip_prefix("P ") {
                    if let Ok(k) = k_str.parse() {
                        m.insert(k, 1);
                    }
                } else if let Some(k_str) = line.strip_prefix("C ") {
                    if let Ok(k) = k_str.parse() {
                        m.compute_if_absent(k, |_| 0);
                    }
                } else if let Some(k_str) = line.strip_prefix("R ") {
                    if let Ok(k) = k_str.parse() {
                        m.remove(&k);
                    }
                }
            }
            if let Some(expected) = expected_keys {
                let actual = m
                    .keys()
                    .map(|k| k.to_string())
                    .collect::<Vec<_>>()
                    .join(" ");
                assert_eq!(actual, expected, "Keys order mismatch for oracle script");
            }
            m
        };

        // Replay ops_1.txt
        let ops_1 = include_str!("../testdata/javamap/ops_1.txt");
        let order_1 = ops_1.lines().last().and_then(|l| l.strip_prefix("ORDER "));
        run_oracle(ops_1, order_1);

        // Replay ops_2.txt
        let ops_2 = include_str!("../testdata/javamap/ops_2.txt");
        let order_2 = ops_2.lines().last().and_then(|l| l.strip_prefix("ORDER "));
        run_oracle(ops_2, order_2);

        // Replay ops_treeify.txt
        let ops_treeify = include_str!("../testdata/javamap/ops_treeify.txt");
        let m = run_oracle(ops_treeify, None);
        assert!(m.treeified(), "ops_treeify must set treeified flag");
    }

    fn replay_oplog_dump_script(script: &str) {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        let mut step = 0i32;
        for line in script.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if line == "NEW" {
                m = JavaHashMap::new();
                continue;
            }
            if let Some(k_str) = line.strip_prefix("P ") {
                let k: i32 = k_str.parse().expect("valid P key");
                m.insert(k, step);
                step += 1;
            } else if let Some(k_str) = line.strip_prefix("G ") {
                let k: i32 = k_str.parse().expect("valid G key");
                let _ = m.get_or_put_with(k, || step);
                step += 1;
            } else if let Some(k_str) = line.strip_prefix("C ") {
                let k: i32 = k_str.parse().expect("valid C key");
                m.compute_if_absent(k, |_| 0);
                step += 1;
            } else if let Some(k_str) = line.strip_prefix("R ") {
                let k: i32 = k_str.parse().expect("valid R key");
                m.remove(&k);
                step += 1;
            } else if let Some(rest) = line.strip_prefix("D ") {
                let (header, order_part) = rest
                    .split_once(" ORDER")
                    .expect("D line must contain ORDER");
                let mut parts = header.split_whitespace();
                let label = parts.next().expect("D label");
                let cap_str = parts
                    .next()
                    .and_then(|s| s.strip_prefix("cap="))
                    .expect("cap=...");
                let expected_cap: usize = cap_str.parse().expect("valid cap");
                let regime = parts
                    .next()
                    .and_then(|s| s.strip_prefix("regime="))
                    .expect("regime=...");
                let expected_order = order_part.trim();
                if regime == "tree" {
                    assert!(
                        m.treeified(),
                        "expected treeified()==true for tree regime on {label}"
                    );
                } else {
                    assert!(
                        !m.treeified(),
                        "expected treeified()==false for list regime on {label}"
                    );
                    assert_eq!(m.capacity(), expected_cap, "capacity mismatch on {label}");
                    let actual_order = m
                        .keys()
                        .map(|k| k.to_string())
                        .collect::<Vec<_>>()
                        .join(" ");
                    assert_eq!(
                        actual_order, expected_order,
                        "key order mismatch on {label}"
                    );
                }
            } else {
                panic!("unknown oplog line: {line}");
            }
        }
    }

    #[test]
    fn op_log_matches_java_oracle_put_threshold() {
        let script = include_str!("../testdata/javamap/ops_put_threshold.txt");
        replay_oplog_dump_script(script);
    }

    #[test]
    fn op_log_matches_java_oracle_getput() {
        let script = include_str!("../testdata/javamap/ops_getput.txt");
        replay_oplog_dump_script(script);
    }

    #[test]
    fn op_log_matches_java_oracle_getput_random() {
        let script = include_str!("../testdata/javamap/ops_getput_random.txt");
        replay_oplog_dump_script(script);
    }

    #[test]
    fn test_javamap_oplog_dump_oracle() {
        replay_oplog_dump_script(include_str!("../testdata/javamap/ops_put_threshold.txt"));
        replay_oplog_dump_script(include_str!("../testdata/javamap/ops_getput.txt"));
        replay_oplog_dump_script(include_str!("../testdata/javamap/ops_getput_random.txt"));
    }

    #[test]
    fn put_treeify_threshold_is_nine_nodes() {
        let mut m16: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..8 {
            m16.insert(i * 16, i);
        }
        assert_eq!(m16.capacity(), 16);
        assert_eq!(
            m16.keys().copied().collect::<Vec<_>>(),
            vec![0, 16, 32, 48, 64, 80, 96, 112]
        );
        m16.insert(8 * 16, 8);
        assert_eq!(m16.capacity(), 32);
        assert_eq!(
            m16.keys().copied().collect::<Vec<_>>(),
            vec![0, 32, 64, 96, 128, 16, 48, 80, 112]
        );

        let mut m32: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..8 {
            m32.insert(i * 32, i);
        }
        assert_eq!(m32.capacity(), 16);
        assert_eq!(
            m32.keys().copied().collect::<Vec<_>>(),
            vec![0, 32, 64, 96, 128, 160, 192, 224]
        );
        m32.insert(8 * 32, 8);
        assert_eq!(m32.capacity(), 32);
        assert_eq!(
            m32.keys().copied().collect::<Vec<_>>(),
            vec![0, 32, 64, 96, 128, 160, 192, 224, 256]
        );
    }

    #[test]
    fn cia_treeify_threshold_is_eight_nodes() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..8 {
            m.compute_if_absent(i * 16, |_| 0);
        }
        assert_eq!(m.capacity(), 32);
        assert_eq!(
            m.keys().copied().collect::<Vec<_>>(),
            vec![96, 64, 32, 0, 112, 80, 48, 16]
        );
    }

    #[test]
    fn get_or_put_matches_java_put_tail_order() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for &k in &[100, 132, 164] {
            m.get_or_put_with(k, || k);
        }
        assert_eq!(m.capacity(), 16);
        assert_eq!(m.keys().copied().collect::<Vec<_>>(), vec![100, 132, 164]);
    }

    #[test]
    fn get_or_put_resizes_after_insert() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for k in 1000..=1012 {
            m.get_or_put_with(k, || k);
        }
        assert_eq!(m.capacity(), 32);
        assert_eq!(
            m.keys().copied().collect::<Vec<_>>(),
            (1000..=1012).collect::<Vec<_>>()
        );
    }

    #[test]
    fn l861q_bin_collision_order() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for k in 55259401..=55259600 {
            m.get_or_put_with(k, || k);
        }
        for k in 55259401..=55259600 {
            if k != 55259524 && k != 55259572 {
                m.remove(&k);
            }
        }
        assert_eq!(
            m.keys().copied().collect::<Vec<_>>(),
            vec![55259524, 55259572]
        );
    }

    #[test]
    fn put_treeified_flag_only_at_nine_nodes_cap64() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        // Grow to capacity 64 using 25 distinct keys in different buckets (threshold 12 -> 24 -> 48)
        for k in 1..=25 {
            m.get_or_put_with(k, || k);
        }
        assert_eq!(m.capacity(), 64);
        assert!(!m.treeified());
        // Insert 8 keys colliding into bucket 0 at capacity 64 (multiples of 64, <= 48 total keys)
        for i in 0..8 {
            m.get_or_put_with(i * 64, || i);
        }
        assert_eq!(m.capacity(), 64);
        assert!(
            !m.treeified(),
            "8 nodes in one bin at cap 64 must not treeify on get_or_put_with"
        );
        m.get_or_put_with(8 * 64, || 8);
        assert!(
            m.treeified(),
            "9th node in one bin at cap 64 must set treeified flag"
        );
    }

    #[test]
    fn get_or_put_matches_java_put_model() {
        fn next_u64(state: &mut u64) -> u64 {
            let mut x = *state;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            *state = x;
            x
        }

        for stream_idx in 0..40u64 {
            let mut rng =
                0xA5A5_5A5A_F0F0_0F0Fu64 ^ (stream_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15));
            if rng == 0 {
                rng = 1;
            }
            let mut arena: JavaHashMap<i32, i32> = JavaHashMap::new();
            let mut legacy: JavaHashMapLegacy<i32, i32> = JavaHashMapLegacy::new();

            for step in 0..2_000usize {
                let r = next_u64(&mut rng);
                let i = ((r >> 8) & 0x3F) as i32;
                let key = match (r >> 16) % 4 {
                    0 => (i << 16) | 0x1234,
                    1 => i * 16,
                    2 => i * 32,
                    _ => i * 65536,
                };
                let val = (step as i32).wrapping_mul(5) ^ key;
                if r % 10 < 7 {
                    let delta = ((r >> 32) & 0xFF) as i32;
                    let a_ref = arena.get_or_put_with(key, || val);
                    *a_ref = a_ref.wrapping_add(delta);
                    let l_ref = legacy.get_or_put_with(key, || val);
                    *l_ref = l_ref.wrapping_add(delta);
                    assert_eq!(*a_ref, *l_ref);
                } else {
                    assert_eq!(arena.remove(&key), legacy.remove(&key));
                }
                assert_eq!(arena.len(), legacy.len());
                assert_eq!(arena.treeified(), legacy.treeified());
                if step % 200 == 0 || step == 1_999 {
                    let a_keys: Vec<i32> = arena.keys().copied().collect();
                    let l_keys: Vec<i32> = legacy.keys().copied().collect();
                    assert_eq!(a_keys, l_keys);
                    let a_into: Vec<(i32, i32)> = arena.clone().into_iter().collect();
                    let l_into: Vec<(i32, i32)> = legacy.clone().into_vec();
                    assert_eq!(a_into, l_into);
                }
            }
        }
    }

    #[test]
    fn op_log_matches_java_oracle_treeify() {
        let ops_treeify = include_str!("../testdata/javamap/ops_treeify.txt");
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for line in ops_treeify.lines() {
            if let Some(k_str) = line.strip_prefix("P ") {
                if let Ok(k) = k_str.parse() {
                    m.insert(k, 1);
                }
            } else if let Some(k_str) = line.strip_prefix("C ") {
                if let Ok(k) = k_str.parse() {
                    m.compute_if_absent(k, |_| 0);
                }
            } else if let Some(k_str) = line.strip_prefix("R ") {
                if let Ok(k) = k_str.parse() {
                    m.remove(&k);
                }
            }
        }
        assert!(m.treeified(), "ops_treeify must set treeified flag");
    }

    #[test]
    fn java_string_hash_golden() {
        assert_eq!(java_string_hash(""), 0);
        assert_eq!(java_string_hash("chr1"), 3_052_836);
        assert_eq!(java_string_hash("chr2"), 3_052_837);
        assert_eq!(java_string_hash("chr3"), 3_052_838);
        assert_eq!(java_string_hash("chr10"), 94_637_964);
        assert_eq!(java_string_hash("chr12"), 94_637_966);
        assert_eq!(java_string_hash("chrX"), 3_052_875);
    }

    #[test]
    fn java_hashmap_key_order_matches_jdk_probe() {
        let k1 = ["chr1", "chrX"];
        assert_eq!(
            java_hashmap_key_order(k1.iter().copied()),
            vec!["chrX", "chr1"]
        );
        let k2 = ["chrX", "chr1"];
        assert_eq!(
            java_hashmap_key_order(k2.iter().copied()),
            vec!["chrX", "chr1"]
        );
        let k3 = ["chr1", "chr12"];
        assert_eq!(
            java_hashmap_key_order(k3.iter().copied()),
            vec!["chr1", "chr12"]
        );
        let k4 = ["chr12", "chr1"];
        assert_eq!(
            java_hashmap_key_order(k4.iter().copied()),
            vec!["chr12", "chr1"]
        );
        // 9 keys colliding into bucket 10 at initial capacity 16 trigger treeifyBin() -> resize()
        // to capacity 32 (where chr1 stays in lo-bucket 10 and chr103 splits into hi-bucket 26).
        let k_collision = [
            "chr103", "chr1", "chr12", "chr23", "chr34", "chr45", "chr56", "chr67", "chr78",
        ];
        let ordered = java_hashmap_key_order(k_collision.iter().copied());
        let pos_chr1 = ordered.iter().position(|&s| s == "chr1");
        let pos_chr103 = ordered.iter().position(|&s| s == "chr103");
        assert!(
            pos_chr1 < pos_chr103,
            "treeifyBin resize from cap 16 to 32 must place chr1 (bucket 10) before chr103 (bucket 26), got {ordered:?}"
        );
        // 11 identical-hash keys ("Aa"/"BB" pairs) force two capacity-doubling resizes
        // (16 -> 32 at key 9, 32 -> 64 at key 10) and then treeify the bucket at key 11,
        // where OpenJDK TreeNode.moveRootToFront moves the red-black root ("AaAaBBBB") to the head.
        let k_treeify = [
            "AaAaAaAa", "AaAaAaBB", "AaAaBBAa", "AaAaBBBB", "AaBBAaAa", "AaBBAaBB", "AaBBBBAa",
            "AaBBBBBB", "BBAaAaAa", "BBAaAaBB", "BBAaBBAa",
        ];
        assert_eq!(
            java_hashmap_key_order(k_treeify.iter().copied()),
            vec![
                "AaAaBBBB", "AaAaAaAa", "AaAaAaBB", "AaAaBBAa", "AaBBAaAa", "AaBBAaBB", "AaBBBBAa",
                "AaBBBBBB", "BBAaAaAa", "BBAaAaBB", "BBAaBBAa",
            ]
        );

        // 12 keys colliding in bin 0 at cap 64 (8 with bit 6 == 0, 4 with bit 6 == 64):
        // - 11th key ("k606") treeifies bin 0;
        // - 12th key ("k727") inserts into the existing TreeNode via putTreeVal and rotates "k727" to root/head;
        // - 37 additional non-bin-0 keys trigger resize 64 -> 128 at size > 48, splitting bin 0 into
        //   retained-tree lo-bin 0 (8 keys -> re-treeified with root "k705" at head) and
        //   untreeified hi-bin 64 (4 keys <= UNTREEIFY_THRESHOLD=6 -> ["k727", "k408", "k529", "k606"]).
        let k_split = [
            "k309", "k507", "k628", "k705", "k749", "k826", "k903", "k947", "k408", "k529", "k606",
            "k727", "pad0", "pad1", "pad2", "pad3", "pad4", "pad5", "pad6", "pad7", "pad8", "pad9",
            "pad10", "pad11", "pad12", "pad13", "pad14", "pad15", "pad16", "pad17", "pad18",
            "pad19", "pad20", "pad21", "pad22", "pad23", "pad24", "pad25", "pad27", "pad28",
            "pad29", "pad30", "pad31", "pad32", "pad33", "pad34", "pad35", "pad36", "pad37",
        ];
        assert_eq!(
            java_hashmap_key_order(k_split.iter().copied()),
            vec![
                "k705", "k309", "k507", "k628", "k749", "k826", "k903", "k947", "pad15", "pad14",
                "pad13", "pad12", "pad11", "pad10", "pad19", "pad18", "pad17", "pad16", "pad7",
                "pad8", "pad9", "pad3", "pad4", "pad5", "pad6", "k727", "k408", "k529", "k606",
                "pad25", "pad24", "pad23", "pad22", "pad21", "pad20", "pad0", "pad1", "pad2",
                "pad29", "pad28", "pad27", "pad37", "pad36", "pad35", "pad34", "pad33", "pad32",
                "pad31", "pad30",
            ]
        );
    }

    #[test]
    fn reserve_nodes_does_not_change_iteration_order() {
        let mut with_reserve: JavaHashMap<i32, i32> = JavaHashMap::new();
        with_reserve.reserve_nodes(512);
        let mut without_reserve: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..512 {
            let k = (i * 37) ^ (i << 7);
            with_reserve.get_or_put_with(k, || i);
            without_reserve.get_or_put_with(k, || i);
        }
        let k1: Vec<(i32, i32)> = with_reserve.into_iter().collect();
        let k2: Vec<(i32, i32)> = without_reserve.into_iter().collect();
        assert_eq!(k1, k2);

        // Opt-in new_key_log captures absent-key allocations across insert / get_or_put_with /
        // compute_if_absent (including free_head reuse and across resize), ignores present-key
        // updates, and does not propagate across Clone.
        let mut logged: JavaHashMap<i32, i32> = JavaHashMap::new();
        assert!(logged.take_new_key_log().is_empty());
        logged.insert(10, 1);
        assert!(logged.take_new_key_log().is_empty());

        logged.enable_new_key_log();
        logged.insert(10, 2); // existing key -> not logged
        logged.insert(20, 3); // new -> logged
        logged.enable_new_key_log(); // idempotent -> preserves prior log entries
        logged.get_or_put_with(20, || 99); // existing -> not logged
        logged.get_or_put_with(30, || 4); // new -> logged
        logged.compute_if_absent(30, |_| 88); // existing -> not logged
        logged.compute_if_absent(40, |_| 5); // new -> logged
        assert_eq!(logged.remove(&20), Some(3));
        logged.get_or_put_with(50, || 6); // reuses free_head -> still logged!
        for k in 100..120 {
            logged.insert(k, k); // triggers resize -> resize must not duplicate log entries
        }

        let mut cloned = logged.clone();
        assert!(cloned.take_new_key_log().is_empty());
        cloned.insert(999, 1);
        assert!(cloned.take_new_key_log().is_empty());

        let mut expected = vec![20, 30, 40, 50];
        expected.extend(100..120);
        assert_eq!(logged.take_new_key_log(), expected);
        // Logging is disabled after take_new_key_log until re-enabled
        logged.insert(777, 1);
        assert!(logged.take_new_key_log().is_empty());
    }

    #[test]
    fn contains_both_and_first_precedes_matches_iter_order() {
        let empty: JavaHashMap<i32, i32> = JavaHashMap::new();
        assert!(!empty.contains_both_and_first_precedes(&1, &1));
        assert!(!empty.contains_both_and_first_precedes(&1, &2));

        let mut map: JavaHashMap<i32, i32> = JavaHashMap::new();
        // Adjacent keys (-1, 0) and (1_048_575, 1_048_576) share bucket 0 under JavaHashMap::spread;
        // 1..500 triggers multiple capacity resizes (16 -> 32 -> 64 -> 128 -> 256 -> 512 -> 1024),
        // and k * 1024 + 7 produces multi-node chains in bucket 7.
        map.insert(-1, -1);
        map.insert(0, 0);
        for k in 1..500 {
            map.insert(k, k * 10);
        }
        for c in 0..7 {
            let colliding = c * 1024 + 7;
            map.insert(colliding, colliding);
        }
        map.compute_if_absent(8 * 1024 + 7, |&k| k);
        assert_eq!(map.remove(&250), Some(2500));
        map.insert(1_048_575, 1);
        map.insert(1_048_576, 2);

        let iter_keys: Vec<i32> = map.iter().map(|(&k, _)| k).collect();
        let rank: std::collections::HashMap<i32, usize> = iter_keys
            .iter()
            .copied()
            .enumerate()
            .map(|(idx, k)| (k, idx))
            .collect();

        for &k in &iter_keys {
            assert!(
                !map.contains_both_and_first_precedes(&k, &k),
                "equal present key {k} must return false"
            );
        }
        for absent in [-999, 250, 500, 999_999, i32::MIN, i32::MAX] {
            assert!(
                !map.contains_both_and_first_precedes(&absent, &absent),
                "equal absent key {absent} must return false"
            );
        }

        for p in (-2..=505).chain(1_048_575..=1_048_577) {
            let prev = p - 1;
            let expected = matches!((rank.get(&prev), rank.get(&p)), (Some(i), Some(j)) if i < j);
            assert_eq!(
                map.contains_both_and_first_precedes(&prev, &p),
                expected,
                "iter relative-order mismatch for ({prev}, {p})"
            );
            let rev_expected =
                matches!((rank.get(&p), rank.get(&prev)), (Some(i), Some(j)) if i < j);
            assert_eq!(
                map.contains_both_and_first_precedes(&p, &prev),
                rev_expected,
                "iter relative-order mismatch for ({p}, {prev})"
            );
        }

        // Verify all pairs across a colliding bucket chain (including absent probe).
        let chain_probes: Vec<i32> = (0..=9).map(|c| c * 1024 + 7).collect();
        for &a in &chain_probes {
            for &b in &chain_probes {
                let expected = matches!((rank.get(&a), rank.get(&b)), (Some(i), Some(j)) if i < j);
                assert_eq!(
                    map.contains_both_and_first_precedes(&a, &b),
                    expected,
                    "chain pair ({a}, {b}) mismatch"
                );
            }
        }

        // Same-bucket adjacent genomic pair (0x000f_ffff, 0x0010_0000) in both insertion orders at cap 16.
        let mut fwd_same_bin: JavaHashMap<i32, i32> = JavaHashMap::new();
        fwd_same_bin.insert(1_048_575, 1);
        fwd_same_bin.insert(1_048_576, 2);
        assert!(fwd_same_bin.contains_both_and_first_precedes(&1_048_575, &1_048_576));
        assert!(!fwd_same_bin.contains_both_and_first_precedes(&1_048_576, &1_048_575));

        let mut rev_same_bin: JavaHashMap<i32, i32> = JavaHashMap::new();
        rev_same_bin.insert(1_048_576, 2);
        rev_same_bin.insert(1_048_575, 1);
        assert!(!rev_same_bin.contains_both_and_first_precedes(&1_048_575, &1_048_576));
        assert!(rev_same_bin.contains_both_and_first_precedes(&1_048_576, &1_048_575));
    }
}
