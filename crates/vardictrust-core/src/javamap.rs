use std::fmt::Debug;
use std::hash::Hash;

/// OpenJDK 8+ `java.util.HashMap` capacity, load-factor, and treeification constants
/// (`HashMap.java:237-275`). Changing any constant breaks `JavaHashMap` iteration parity.
const DEFAULT_CAPACITY: usize = 16;
const MAXIMUM_CAPACITY: usize = 1 << 30;
const TREEIFY_THRESHOLD: usize = 8;
const UNTREEIFY_THRESHOLD: usize = 6;
const MIN_TREEIFY_CAPACITY: usize = 64;
const DEFAULT_LOAD_FACTOR: f32 = 0.75;
const NONE: u32 = u32::MAX;

#[inline]
fn next_resize_plan(old_cap: usize) -> Option<(usize, usize)> {
    if old_cap >= MAXIMUM_CAPACITY {
        return None;
    }
    let new_cap = (old_cap << 1).min(MAXIMUM_CAPACITY);
    let new_thr = if new_cap >= MAXIMUM_CAPACITY {
        i32::MAX as usize
    } else {
        let ft = (new_cap as f32) * DEFAULT_LOAD_FACTOR;
        if ft < (MAXIMUM_CAPACITY as f32) {
            ft as usize
        } else {
            i32::MAX as usize
        }
    };
    Some((new_cap, new_thr))
}

#[derive(Clone, Debug)]
struct Node<K, V> {
    key: K,
    val: V,
    next: u32,
}

/// Cold sidecar Red-Black tree + doubly-linked list links corresponding to
/// `java.util.HashMap.TreeNode` (`HashMap.java:1967-1975`). Allocated lazily only when
/// at least one bucket crosses `TREEIFY_THRESHOLD` at `capacity >= MIN_TREEIFY_CAPACITY`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct TreeLink {
    parent: u32,
    left: u32,
    right: u32,
    prev: u32,
    red: bool,
    is_tree: bool,
}

impl Default for TreeLink {
    #[inline]
    fn default() -> Self {
        Self {
            parent: NONE,
            left: NONE,
            right: NONE,
            prev: NONE,
            red: false,
            is_tree: false,
        }
    }
}

enum TreeProbe {
    Found(u32),
    Vacant { root: u32, xp: u32, dir: i32 },
}

#[derive(Debug)]
#[expect(
    clippy::box_collection,
    reason = "8-byte niche Option<Box<Vec<K>>> keeps JavaHashMap at 96B when key logging is off"
)]
struct NewKeyLog<K>(Option<Box<Vec<K>>>);

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

#[cold]
#[inline(never)]
fn fallback_slot_mut<V: Default>(slot: &mut Option<Box<V>>) -> &mut V {
    slot.get_or_insert_with(|| Box::new(V::default()))
}

#[derive(Clone, Debug)]
pub struct JavaHashMap<K, V> {
    buckets: Vec<u32>,
    nodes: Vec<Node<K, V>>,
    #[expect(
        clippy::box_collection,
        reason = "8-byte niche Option<Box<Vec<TreeLink>>> avoids +16B per JavaHashMap on cold non-treeifying path"
    )]
    tree_links: Option<Box<Vec<TreeLink>>>,
    free_head: u32,
    active_tree_bins: u32,
    size: usize,
    threshold: usize,
    fallback: Option<Box<V>>,
    new_key_log: NewKeyLog<K>,
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
#[derive(Debug)]
pub struct Iter<'a, K, V> {
    buckets: &'a [u32],
    nodes: &'a [Node<K, V>],
    bucket_idx: usize,
    curr: u32,
}

impl<K, V> Clone for Iter<'_, K, V> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            buckets: self.buckets,
            nodes: self.nodes,
            bucket_idx: self.bucket_idx,
            curr: self.curr,
        }
    }
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

mod private {
    pub trait Sealed {}
    impl Sealed for i32 {}
    impl Sealed for str {}
    impl Sealed for String {}
    impl<T: Sealed + ?Sized> Sealed for &T {}
}

/// Computes the 32-bit signed `Object.hashCode()` and `Comparable.compareTo()` total order
/// used by [`JavaHashMap`]. Sealed to Java `Comparable` key types (`i32`, `str`, `String`, `&T`).
pub trait JavaHash: private::Sealed {
    fn java_hash(&self) -> i32;
    fn java_compare(&self, other: &Self) -> std::cmp::Ordering;
}

impl JavaHash for i32 {
    #[inline]
    fn java_hash(&self) -> i32 {
        *self
    }

    #[inline]
    fn java_compare(&self, other: &Self) -> std::cmp::Ordering {
        self.cmp(other)
    }
}

impl JavaHash for str {
    #[inline]
    fn java_hash(&self) -> i32 {
        java_string_hash(self)
    }

    #[inline]
    fn java_compare(&self, other: &Self) -> std::cmp::Ordering {
        if self.is_ascii() && other.is_ascii() {
            return self.as_bytes().cmp(other.as_bytes());
        }
        self.encode_utf16().cmp(other.encode_utf16())
    }
}

impl JavaHash for String {
    #[inline]
    fn java_hash(&self) -> i32 {
        java_string_hash(self)
    }

    #[inline]
    fn java_compare(&self, other: &Self) -> std::cmp::Ordering {
        self.as_str().java_compare(other.as_str())
    }
}

impl<T: JavaHash + ?Sized> JavaHash for &T {
    #[inline]
    fn java_hash(&self) -> i32 {
        (**self).java_hash()
    }

    #[inline]
    fn java_compare(&self, other: &Self) -> std::cmp::Ordering {
        (**self).java_compare(*other)
    }
}

impl<K, V> JavaHashMap<K, V>
where
    K: Clone + JavaHash + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    #[must_use]
    pub fn new() -> Self {
        Self {
            buckets: Vec::new(),
            nodes: Vec::new(),
            tree_links: None,
            free_head: NONE,
            active_tree_bins: 0,
            size: 0,
            threshold: 0,
            fallback: None,
            new_key_log: NewKeyLog::default(),
        }
    }

    /// Records newly allocated keys during an insert-only window (`sv::partial_pipeline`).
    /// Idempotent when already enabled; `remove` within the window is outside the contract.
    pub(crate) fn enable_new_key_log(&mut self) {
        self.new_key_log
            .0
            .get_or_insert_with(|| Box::new(Vec::new()));
    }

    pub(crate) fn take_new_key_log(&mut self) -> Vec<K> {
        self.new_key_log.0.take().map(|b| *b).unwrap_or_default()
    }

    #[inline]
    fn spread(key: &K) -> u32 {
        let h = key.java_hash() as u32;
        h ^ (h >> 16)
    }

    /// Signed 32-bit `TreeNode` comparison (`HashMap.java:2023-2035, 2089-2096, 2141-2160`):
    /// compares `ph` vs `h` as signed `i32`, breaking hash ties via `Comparable.compareTo`
    /// (`JavaHash::java_compare`).
    #[inline]
    fn tree_cmp(h: i32, key: &K, ph: i32, pk: &K) -> i32 {
        if ph > h {
            -1
        } else if ph < h {
            1
        } else {
            match key.java_compare(pk) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            }
        }
    }

    #[inline]
    fn is_tree_node(&self, idx: u32) -> bool {
        self.tree_links
            .as_deref()
            .and_then(|links| links.get(idx as usize))
            .is_some_and(|l| l.is_tree)
    }

    #[inline]
    fn link(&self, idx: u32) -> Option<&TreeLink> {
        self.tree_links
            .as_deref()
            .and_then(|links| links.get(idx as usize))
    }

    #[inline]
    fn link_mut(&mut self, idx: u32) -> Option<&mut TreeLink> {
        self.tree_links
            .as_deref_mut()
            .and_then(|links| links.get_mut(idx as usize))
    }

    #[inline]
    fn reset_link(&mut self, idx: u32) {
        if let Some(slot) = self.link_mut(idx) {
            *slot = TreeLink::default();
        }
    }

    fn ensure_tree_links(&mut self) {
        let target_len = self.nodes.len();
        match self.tree_links.as_deref_mut() {
            Some(links) => {
                if links.len() < target_len {
                    links.resize(target_len, TreeLink::default());
                }
            }
            None => {
                self.tree_links = Some(Box::new(vec![TreeLink::default(); target_len]));
            }
        }
    }

    fn alloc_node(&mut self, key: K, val: V, next: u32) -> u32 {
        if self.free_head == NONE && self.nodes.len() >= (u32::MAX as usize) {
            return NONE;
        }
        if let Some(log) = self.new_key_log.0.as_deref_mut() {
            log.push(key.clone());
        }
        if self.free_head != NONE {
            let idx = self.free_head;
            if let Some(slot) = self.nodes.get_mut(idx as usize) {
                self.free_head = slot.next;
                *slot = Node { key, val, next };
                self.reset_link(idx);
                return idx;
            }
            self.free_head = NONE;
            if self.nodes.len() >= (u32::MAX as usize) {
                return NONE;
            }
        }
        let idx = self.nodes.len() as u32;
        self.nodes.push(Node { key, val, next });
        if let Some(links) = self.tree_links.as_deref_mut() {
            links.push(TreeLink::default());
        }
        idx
    }

    #[must_use]
    pub fn treeified(&self) -> bool {
        self.active_tree_bins > 0
    }

    #[must_use]
    pub fn capacity(&self) -> usize {
        self.buckets.len()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.size
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// Pre-allocates capacity in the backing `nodes` arena without altering bucket
    /// count, resize threshold, or Java `HashMap` bucket-chain iteration order (`E4`).
    pub fn reserve_nodes(&mut self, additional: usize) {
        self.nodes.reserve(additional);
        if let Some(links) = self.tree_links.as_deref_mut() {
            links.reserve(additional);
        }
    }

    #[must_use]
    pub fn contains_key(&self, key: &K) -> bool {
        self.find_node(key) != NONE
    }

    /// Returns `true` iff both `first` and `second` are present in the map and `first`
    /// strictly precedes `second` in [`Self::iter`] order (`false` when `first == second`).
    #[must_use]
    pub fn contains_both_and_first_precedes(&self, first: &K, second: &K) -> bool {
        if first == second || self.buckets.is_empty() {
            return false;
        }
        let n = self.buckets.len();
        let h_first = Self::spread(first);
        let h_second = Self::spread(second);
        let idx_first = (h_first as usize) & (n - 1);
        let idx_second = (h_second as usize) & (n - 1);
        match idx_first.cmp(&idx_second) {
            std::cmp::Ordering::Greater => false,
            std::cmp::Ordering::Less => self.contains_key(first) && self.contains_key(second),
            std::cmp::Ordering::Equal => {
                let head = self.buckets.get(idx_first).copied().unwrap_or(NONE);
                if head == NONE {
                    return false;
                }
                if self.is_tree_node(head) {
                    if self.get_tree_node(head, h_first as i32, first) == NONE
                        || self.get_tree_node(head, h_second as i32, second) == NONE
                    {
                        return false;
                    }
                    let mut curr = head;
                    for _ in 0..=self.nodes.len() {
                        if curr == NONE {
                            break;
                        }
                        let Some(node) = self.nodes.get(curr as usize) else {
                            return false;
                        };
                        if &node.key == first {
                            return true;
                        }
                        if &node.key == second {
                            return false;
                        }
                        curr = node.next;
                    }
                    false
                } else {
                    let mut curr = head;
                    let mut saw_first = false;
                    for _ in 0..=self.nodes.len() {
                        if curr == NONE {
                            break;
                        }
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
    }

    fn find_node(&self, key: &K) -> u32 {
        if self.buckets.is_empty() {
            return NONE;
        }
        let hash = Self::spread(key);
        let n = self.buckets.len();
        let idx = (hash as usize) & (n - 1);
        let first = self.buckets.get(idx).copied().unwrap_or(NONE);
        if first == NONE {
            return NONE;
        }
        let Some(first_node) = self.nodes.get(first as usize) else {
            return NONE;
        };
        if &first_node.key == key {
            return first;
        }
        let mut curr = first_node.next;
        if curr == NONE {
            return NONE;
        }
        if self.is_tree_node(first) {
            return self.get_tree_node(first, hash as i32, key);
        }
        for _ in 0..=self.nodes.len() {
            if curr == NONE {
                break;
            }
            let Some(node) = self.nodes.get(curr as usize) else {
                return NONE;
            };
            if &node.key == key {
                return curr;
            }
            curr = node.next;
        }
        NONE
    }

    #[must_use]
    pub fn get(&self, key: &K) -> Option<&V> {
        let idx = self.find_node(key);
        if idx == NONE {
            None
        } else {
            self.nodes.get(idx as usize).map(|n| &n.val)
        }
    }

    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        let idx = self.find_node(key);
        if idx == NONE {
            None
        } else {
            self.nodes.get_mut(idx as usize).map(|n| &mut n.val)
        }
    }

    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if self.buckets.is_empty() {
            self.resize();
        }
        let hash = Self::spread(&key);
        let n = self.buckets.len();
        if n == 0 {
            return None;
        }
        let idx = (hash as usize) & (n - 1);
        let first = self.buckets.get(idx).copied().unwrap_or(NONE);

        if first == NONE {
            let new_idx = self.alloc_node(key, value, NONE);
            if new_idx == NONE {
                return None;
            }
            if let Some(b) = self.buckets.get_mut(idx) {
                *b = new_idx;
            }
        } else if self
            .nodes
            .get(first as usize)
            .is_some_and(|nd| nd.key == key)
        {
            if let Some(node_mut) = self.nodes.get_mut(first as usize) {
                return Some(std::mem::replace(&mut node_mut.val, value));
            }
            return None;
        } else if self.is_tree_node(first) {
            match self.probe_tree_slot(first, hash as i32, &key) {
                TreeProbe::Found(e) => {
                    if let Some(node_mut) = self.nodes.get_mut(e as usize) {
                        return Some(std::mem::replace(&mut node_mut.val, value));
                    }
                    return None;
                }
                TreeProbe::Vacant { root, xp, dir } => {
                    if self.insert_tree_leaf(idx, root, xp, dir, key, value) == NONE {
                        return None;
                    }
                }
            }
        } else {
            let mut p = first;
            for bin_count in 0..=self.nodes.len() {
                let e = self.nodes.get(p as usize).map_or(NONE, |nd| nd.next);
                if e == NONE {
                    let new_idx = self.alloc_node(key, value, NONE);
                    if new_idx == NONE {
                        return None;
                    }
                    if let Some(p_mut) = self.nodes.get_mut(p as usize) {
                        p_mut.next = new_idx;
                    }
                    if bin_count >= TREEIFY_THRESHOLD - 1 {
                        self.treeify_bin(hash);
                    }
                    break;
                }
                if self.nodes.get(e as usize).is_some_and(|nd| nd.key == key) {
                    if let Some(node_mut) = self.nodes.get_mut(e as usize) {
                        return Some(std::mem::replace(&mut node_mut.val, value));
                    }
                    return None;
                }
                p = e;
            }
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
        let hash = Self::spread(&key);
        let n = self.buckets.len();
        if n == 0 {
            return fallback_slot_mut(&mut self.fallback);
        }
        let idx = (hash as usize) & (n - 1);
        let first = self.buckets.get(idx).copied().unwrap_or(NONE);

        let new_idx = if first == NONE {
            let ni = self.alloc_node(key, f(), NONE);
            if ni == NONE {
                return fallback_slot_mut(&mut self.fallback);
            }
            if let Some(b) = self.buckets.get_mut(idx) {
                *b = ni;
            }
            ni
        } else if self
            .nodes
            .get(first as usize)
            .is_some_and(|nd| nd.key == key)
        {
            return self
                .nodes
                .get_mut(first as usize)
                .map(|nd| &mut nd.val)
                .unwrap_or_else(|| fallback_slot_mut(&mut self.fallback));
        } else if self.is_tree_node(first) {
            match self.probe_tree_slot(first, hash as i32, &key) {
                TreeProbe::Found(e) => {
                    return self
                        .nodes
                        .get_mut(e as usize)
                        .map(|nd| &mut nd.val)
                        .unwrap_or_else(|| fallback_slot_mut(&mut self.fallback));
                }
                TreeProbe::Vacant { root, xp, dir } => {
                    let ni = self.insert_tree_leaf(idx, root, xp, dir, key, f());
                    if ni == NONE {
                        return fallback_slot_mut(&mut self.fallback);
                    }
                    ni
                }
            }
        } else {
            let mut p = first;
            let mut inserted = NONE;
            for bin_count in 0..=self.nodes.len() {
                let e = self.nodes.get(p as usize).map_or(NONE, |nd| nd.next);
                if e == NONE {
                    let ni = self.alloc_node(key, f(), NONE);
                    if ni == NONE {
                        return fallback_slot_mut(&mut self.fallback);
                    }
                    if let Some(p_mut) = self.nodes.get_mut(p as usize) {
                        p_mut.next = ni;
                    }
                    if bin_count >= TREEIFY_THRESHOLD - 1 {
                        self.treeify_bin(hash);
                    }
                    inserted = ni;
                    break;
                }
                if self.nodes.get(e as usize).is_some_and(|nd| nd.key == key) {
                    return self
                        .nodes
                        .get_mut(e as usize)
                        .map(|nd| &mut nd.val)
                        .unwrap_or_else(|| fallback_slot_mut(&mut self.fallback));
                }
                p = e;
            }
            if inserted == NONE {
                return fallback_slot_mut(&mut self.fallback);
            }
            inserted
        };

        self.size += 1;
        if self.size > self.threshold {
            self.resize();
        }

        if let Some(nd) = self.nodes.get_mut(new_idx as usize) {
            &mut nd.val
        } else {
            fallback_slot_mut(&mut self.fallback)
        }
    }

    /// Emulates OpenJDK `HashMap.computeIfAbsent` (`HashMap.java:1196-1249`):
    /// pre-lookup resize (`size > threshold`), `TreeNode.getTreeNode`/`putTreeVal` when
    /// bucket is treeified, head-prepend + `treeifyBin` at `binCount >= 7` (8th node)
    /// when plain list, and no post-insert resize.
    #[deprecated(
        since = "1.0.10",
        note = "use JavaHashMap::get_or_put_with for get+put semantics; compute_if_absent emulates OpenJDK HashMap.computeIfAbsent"
    )]
    pub fn compute_if_absent<F>(&mut self, key: K, f: F) -> &mut V
    where
        F: FnOnce(&K) -> V,
    {
        if self.size > self.threshold || self.buckets.is_empty() {
            self.resize();
        }
        let hash = Self::spread(&key);
        let n = self.buckets.len();
        if n == 0 {
            return fallback_slot_mut(&mut self.fallback);
        }
        let idx = (hash as usize) & (n - 1);
        let first = self.buckets.get(idx).copied().unwrap_or(NONE);

        let mut bin_count = 0usize;
        let mut tree_vacant: Option<(u32, u32, i32)> = None;
        if first != NONE {
            if self.is_tree_node(first) {
                match self.probe_tree_slot(first, hash as i32, &key) {
                    TreeProbe::Found(old) => {
                        return self
                            .nodes
                            .get_mut(old as usize)
                            .map(|nd| &mut nd.val)
                            .unwrap_or_else(|| fallback_slot_mut(&mut self.fallback));
                    }
                    TreeProbe::Vacant { root, xp, dir } => {
                        tree_vacant = Some((root, xp, dir));
                    }
                }
            } else {
                let mut e = first;
                let mut old = NONE;
                for _ in 0..=self.nodes.len() {
                    if e == NONE {
                        break;
                    }
                    let Some(node) = self.nodes.get(e as usize) else {
                        break;
                    };
                    if node.key == key {
                        old = e;
                        break;
                    }
                    bin_count += 1;
                    e = node.next;
                }
                if old != NONE {
                    return self
                        .nodes
                        .get_mut(old as usize)
                        .map(|nd| &mut nd.val)
                        .unwrap_or_else(|| fallback_slot_mut(&mut self.fallback));
                }
            }
        }

        let val = f(&key);
        let target_idx = if let Some((root, xp, dir)) = tree_vacant {
            let ni = self.insert_tree_leaf(idx, root, xp, dir, key, val);
            if ni != NONE {
                self.size += 1;
            }
            ni
        } else {
            let new_idx = self.alloc_node(key, val, first);
            if new_idx != NONE {
                if let Some(b) = self.buckets.get_mut(idx) {
                    *b = new_idx;
                }
                if bin_count >= TREEIFY_THRESHOLD - 1 {
                    self.treeify_bin(hash);
                }
                self.size += 1;
            }
            new_idx
        };

        if let Some(nd) = self.nodes.get_mut(target_idx as usize) {
            &mut nd.val
        } else {
            fallback_slot_mut(&mut self.fallback)
        }
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
        if self.buckets.is_empty() {
            return None;
        }
        let hash = Self::spread(key);
        let n = self.buckets.len();
        let idx = (hash as usize) & (n - 1);
        let p = self.buckets.get(idx).copied().unwrap_or(NONE);
        if p == NONE {
            return None;
        }

        let (node, prev_list) = if self.nodes.get(p as usize).is_some_and(|nd| &nd.key == key) {
            (p, NONE)
        } else {
            let e0 = self.nodes.get(p as usize).map_or(NONE, |nd| nd.next);
            if e0 == NONE {
                return None;
            }
            if self.is_tree_node(p) {
                let found = self.get_tree_node(p, hash as i32, key);
                if found == NONE {
                    return None;
                }
                (found, NONE)
            } else {
                let mut prev = p;
                let mut e = e0;
                let mut found = NONE;
                for _ in 0..=self.nodes.len() {
                    if e == NONE {
                        break;
                    }
                    let Some(nd) = self.nodes.get(e as usize) else {
                        break;
                    };
                    if &nd.key == key {
                        found = e;
                        break;
                    }
                    prev = e;
                    e = nd.next;
                }
                if found == NONE {
                    return None;
                }
                (found, prev)
            }
        };

        if self.is_tree_node(node) {
            self.remove_tree_node(idx, node, true);
        } else if node == p {
            let next = self.nodes.get(node as usize).map_or(NONE, |nd| nd.next);
            if let Some(b) = self.buckets.get_mut(idx) {
                *b = next;
            }
        } else {
            let next = self.nodes.get(node as usize).map_or(NONE, |nd| nd.next);
            if let Some(pn) = self.nodes.get_mut(prev_list as usize) {
                pn.next = next;
            }
        }

        self.size = self.size.saturating_sub(1);
        self.reset_link(node);
        if let Some(removed_node) = self.nodes.get_mut(node as usize) {
            removed_node.next = self.free_head;
            self.free_head = node;
            return Some(std::mem::take(&mut removed_node.val));
        }
        None
    }

    /// Ports `HashMap.treeifyBin` (`HashMap.java:762-781`).
    fn treeify_bin(&mut self, hash: u32) {
        let n = self.buckets.len();
        if n < MIN_TREEIFY_CAPACITY {
            self.resize();
            return;
        }
        let index = (hash as usize) & (n - 1);
        let head = self.buckets.get(index).copied().unwrap_or(NONE);
        if head != NONE && !self.is_tree_node(head) {
            self.ensure_tree_links();
            self.treeify_bucket(index);
            self.active_tree_bins = self.active_tree_bins.saturating_add(1);
        }
    }

    /// Ports `TreeNode.root()` (`HashMap.java:1980-1986`).
    fn tree_root(&self, start: u32) -> u32 {
        let Some(links) = self.tree_links.as_deref() else {
            return start;
        };
        let mut r = start;
        for _ in 0..=links.len() {
            let Some(link) = links.get(r as usize) else {
                return r;
            };
            if link.parent == NONE {
                return r;
            }
            r = link.parent;
        }
        r
    }

    /// Ports `TreeNode.moveRootToFront` (`HashMap.java:1991-2011`).
    fn move_root_to_front(&mut self, bucket_idx: usize, root: u32) {
        let first = self.buckets.get(bucket_idx).copied().unwrap_or(NONE);
        if root != NONE && root != first {
            if let Some(b) = self.buckets.get_mut(bucket_idx) {
                *b = root;
            }
            let rp = self.link(root).map_or(NONE, |l| l.prev);
            let rn = self.nodes.get(root as usize).map_or(NONE, |n| n.next);
            if rn != NONE {
                if let Some(rnl) = self.link_mut(rn) {
                    rnl.prev = rp;
                }
            }
            if rp != NONE {
                if let Some(rpn) = self.nodes.get_mut(rp as usize) {
                    rpn.next = rn;
                }
            }
            if first != NONE {
                if let Some(fl) = self.link_mut(first) {
                    fl.prev = root;
                }
            }
            if let Some(rn_mut) = self.nodes.get_mut(root as usize) {
                rn_mut.next = first;
            }
            if let Some(rl_mut) = self.link_mut(root) {
                rl_mut.prev = NONE;
            }
        }
    }

    /// Ports `TreeNode.getTreeNode` + `TreeNode.find` (`HashMap.java:2018-2050`).
    fn get_tree_node(&self, first: u32, h: i32, key: &K) -> u32 {
        let root = if self.link(first).is_some_and(|l| l.parent != NONE) {
            self.tree_root(first)
        } else {
            first
        };
        let Some(links) = self.tree_links.as_deref() else {
            return NONE;
        };
        let mut p = root;
        for _ in 0..=self.nodes.len() {
            if p == NONE {
                return NONE;
            }
            let (Some(node), Some(link)) = (self.nodes.get(p as usize), links.get(p as usize))
            else {
                return NONE;
            };
            let ph = Self::spread(&node.key) as i32;
            let pl = link.left;
            let pr = link.right;
            if ph > h {
                p = pl;
            } else if ph < h {
                p = pr;
            } else if &node.key == key {
                return p;
            } else if pl == NONE {
                p = pr;
            } else if pr == NONE {
                p = pl;
            } else {
                p = match key.java_compare(&node.key) {
                    std::cmp::Ordering::Less => pl,
                    std::cmp::Ordering::Equal => return p,
                    std::cmp::Ordering::Greater => pr,
                };
            }
        }
        NONE
    }

    /// Ports the search phase of `TreeNode.putTreeVal` (`HashMap.java:2134-2163`).
    fn probe_tree_slot(&self, first: u32, h: i32, key: &K) -> TreeProbe {
        let root = if self.link(first).is_some_and(|l| l.parent != NONE) {
            self.tree_root(first)
        } else {
            first
        };
        let Some(links) = self.tree_links.as_deref() else {
            return TreeProbe::Vacant {
                root,
                xp: first,
                dir: 1,
            };
        };
        let mut p = root;
        for _ in 0..=self.nodes.len() {
            let (Some(node), Some(link)) = (self.nodes.get(p as usize), links.get(p as usize))
            else {
                break;
            };
            let ph = Self::spread(&node.key) as i32;
            if ph == h && &node.key == key {
                return TreeProbe::Found(p);
            }
            let dir = Self::tree_cmp(h, key, ph, &node.key);
            if dir == 0 {
                return TreeProbe::Found(p);
            }
            let xp = p;
            p = if dir <= 0 { link.left } else { link.right };
            if p == NONE {
                return TreeProbe::Vacant { root, xp, dir };
            }
        }
        TreeProbe::Vacant {
            root,
            xp: first,
            dir: 1,
        }
    }

    /// Ports the leaf-link, balance, and `moveRootToFront` phase of `TreeNode.putTreeVal`
    /// (`HashMap.java:2163-2176`).
    fn insert_tree_leaf(
        &mut self,
        bucket_idx: usize,
        root: u32,
        xp: u32,
        dir: i32,
        key: K,
        val: V,
    ) -> u32 {
        let xpn = self.nodes.get(xp as usize).map_or(NONE, |nd| nd.next);
        let x = self.alloc_node(key, val, xpn);
        if x == NONE {
            return NONE;
        }
        self.ensure_tree_links();
        if let Some(xp_link) = self.link_mut(xp) {
            if dir <= 0 {
                xp_link.left = x;
            } else {
                xp_link.right = x;
            }
        }
        if let Some(xp_node) = self.nodes.get_mut(xp as usize) {
            xp_node.next = x;
        }
        if let Some(x_link) = self.link_mut(x) {
            *x_link = TreeLink {
                parent: xp,
                left: NONE,
                right: NONE,
                prev: xp,
                red: true,
                is_tree: true,
            };
        }
        if xpn != NONE {
            if let Some(xpn_link) = self.link_mut(xpn) {
                xpn_link.prev = x;
            }
        }
        let new_root = self.balance_insertion(root, x);
        self.move_root_to_front(bucket_idx, new_root);
        x
    }

    /// Ports `TreeNode.treeify` (`HashMap.java:2072-2112`) on bucket `bucket_idx`.
    /// Caller updates `active_tree_bins`.
    fn treeify_bucket(&mut self, bucket_idx: usize) {
        let head = self.buckets.get(bucket_idx).copied().unwrap_or(NONE);
        if head == NONE {
            return;
        }
        self.ensure_tree_links();
        let max_steps = self.nodes.len();

        let mut tl = NONE;
        let mut curr = head;
        for _ in 0..=max_steps {
            if curr == NONE {
                break;
            }
            let next = self.nodes.get(curr as usize).map_or(NONE, |nd| nd.next);
            if let Some(link) = self.link_mut(curr) {
                *link = TreeLink {
                    parent: NONE,
                    left: NONE,
                    right: NONE,
                    prev: tl,
                    red: false,
                    is_tree: true,
                };
            }
            tl = curr;
            curr = next;
        }

        let mut root = NONE;
        let mut x = head;
        for _ in 0..=max_steps {
            if x == NONE {
                break;
            }
            let next = self.nodes.get(x as usize).map_or(NONE, |nd| nd.next);
            if let Some(xl) = self.link_mut(x) {
                xl.left = NONE;
                xl.right = NONE;
            }
            if root == NONE {
                if let Some(xl) = self.link_mut(x) {
                    xl.parent = NONE;
                    xl.red = false;
                }
                root = x;
            } else {
                let mut p = root;
                for _ in 0..=max_steps {
                    let (dir, pl, pr) = {
                        let (Some(xn), Some(pn), Some(plk)) = (
                            self.nodes.get(x as usize),
                            self.nodes.get(p as usize),
                            self.link(p),
                        ) else {
                            break;
                        };
                        let h = Self::spread(&xn.key) as i32;
                        let ph = Self::spread(&pn.key) as i32;
                        (Self::tree_cmp(h, &xn.key, ph, &pn.key), plk.left, plk.right)
                    };
                    let xp = p;
                    p = if dir <= 0 { pl } else { pr };
                    if p == NONE {
                        if let Some(xl) = self.link_mut(x) {
                            xl.parent = xp;
                        }
                        if let Some(xpl) = self.link_mut(xp) {
                            if dir <= 0 {
                                xpl.left = x;
                            } else {
                                xpl.right = x;
                            }
                        }
                        root = self.balance_insertion(root, x);
                        break;
                    }
                }
            }
            x = next;
        }
        self.move_root_to_front(bucket_idx, root);
    }

    /// Ports `TreeNode.untreeify` (`HashMap.java:2118-2129`), resetting all `TreeLink`
    /// entries in the chain starting at `head` back to `TreeLink::default()`.
    fn untreeify_chain(&mut self, head: u32) {
        let mut q = head;
        for _ in 0..=self.nodes.len() {
            if q == NONE {
                break;
            }
            let next = self.nodes.get(q as usize).map_or(NONE, |nd| nd.next);
            self.reset_link(q);
            q = next;
        }
    }

    /// Ports `TreeNode.rotateLeft` (`HashMap.java:2348-2364`).
    fn rotate_left(&mut self, mut root: u32, p: u32) -> u32 {
        if p == NONE {
            return root;
        }
        let r = self.link(p).map_or(NONE, |l| l.right);
        if r == NONE {
            return root;
        }
        let rl = self.link(r).map_or(NONE, |l| l.left);
        if let Some(pl) = self.link_mut(p) {
            pl.right = rl;
        }
        if rl != NONE {
            if let Some(rll) = self.link_mut(rl) {
                rll.parent = p;
            }
        }
        let pp = self.link(p).map_or(NONE, |l| l.parent);
        if let Some(rlk) = self.link_mut(r) {
            rlk.parent = pp;
        }
        if pp == NONE {
            root = r;
            if let Some(rlk) = self.link_mut(r) {
                rlk.red = false;
            }
        } else if self.link(pp).is_some_and(|l| l.left == p) {
            if let Some(ppl) = self.link_mut(pp) {
                ppl.left = r;
            }
        } else if let Some(ppl) = self.link_mut(pp) {
            ppl.right = r;
        }
        if let Some(rlk) = self.link_mut(r) {
            rlk.left = p;
        }
        if let Some(pl) = self.link_mut(p) {
            pl.parent = r;
        }
        root
    }

    /// Ports `TreeNode.rotateRight` (`HashMap.java:2366-2382`).
    fn rotate_right(&mut self, mut root: u32, p: u32) -> u32 {
        if p == NONE {
            return root;
        }
        let l = self.link(p).map_or(NONE, |lk| lk.left);
        if l == NONE {
            return root;
        }
        let lr = self.link(l).map_or(NONE, |lk| lk.right);
        if let Some(pl) = self.link_mut(p) {
            pl.left = lr;
        }
        if lr != NONE {
            if let Some(lrl) = self.link_mut(lr) {
                lrl.parent = p;
            }
        }
        let pp = self.link(p).map_or(NONE, |lk| lk.parent);
        if let Some(llk) = self.link_mut(l) {
            llk.parent = pp;
        }
        if pp == NONE {
            root = l;
            if let Some(llk) = self.link_mut(l) {
                llk.red = false;
            }
        } else if self.link(pp).is_some_and(|lk| lk.right == p) {
            if let Some(ppl) = self.link_mut(pp) {
                ppl.right = l;
            }
        } else if let Some(ppl) = self.link_mut(pp) {
            ppl.left = l;
        }
        if let Some(llk) = self.link_mut(l) {
            llk.right = p;
        }
        if let Some(pl) = self.link_mut(p) {
            pl.parent = l;
        }
        root
    }

    /// Ports `TreeNode.balanceInsertion` (`HashMap.java:2384-2437`).
    fn balance_insertion(&mut self, mut root: u32, mut x: u32) -> u32 {
        if let Some(xl) = self.link_mut(x) {
            xl.red = true;
        }
        for _ in 0..=self.nodes.len() {
            let xp = self.link(x).map_or(NONE, |l| l.parent);
            if xp == NONE {
                if let Some(xl) = self.link_mut(x) {
                    xl.red = false;
                }
                return x;
            }
            let (xp_red, xpp) = self.link(xp).map_or((false, NONE), |l| (l.red, l.parent));
            if !xp_red || xpp == NONE {
                return root;
            }
            let (xppl, xppr) = self.link(xpp).map_or((NONE, NONE), |l| (l.left, l.right));
            if xp == xppl {
                if xppr != NONE && self.link(xppr).is_some_and(|l| l.red) {
                    if let Some(l) = self.link_mut(xppr) {
                        l.red = false;
                    }
                    if let Some(l) = self.link_mut(xp) {
                        l.red = false;
                    }
                    if let Some(l) = self.link_mut(xpp) {
                        l.red = true;
                    }
                    x = xpp;
                } else {
                    if self.link(xp).is_some_and(|l| x == l.right) {
                        x = xp;
                        root = self.rotate_left(root, x);
                    }
                    let xp = self.link(x).map_or(NONE, |l| l.parent);
                    let xpp = self.link(xp).map_or(NONE, |l| l.parent);
                    if xp != NONE {
                        if let Some(l) = self.link_mut(xp) {
                            l.red = false;
                        }
                        if xpp != NONE {
                            if let Some(l) = self.link_mut(xpp) {
                                l.red = true;
                            }
                            root = self.rotate_right(root, xpp);
                        }
                    }
                }
            } else if xppl != NONE && self.link(xppl).is_some_and(|l| l.red) {
                if let Some(l) = self.link_mut(xppl) {
                    l.red = false;
                }
                if let Some(l) = self.link_mut(xp) {
                    l.red = false;
                }
                if let Some(l) = self.link_mut(xpp) {
                    l.red = true;
                }
                x = xpp;
            } else {
                if self.link(xp).is_some_and(|l| x == l.left) {
                    x = xp;
                    root = self.rotate_right(root, x);
                }
                let xp = self.link(x).map_or(NONE, |l| l.parent);
                let xpp = self.link(xp).map_or(NONE, |l| l.parent);
                if xp != NONE {
                    if let Some(l) = self.link_mut(xp) {
                        l.red = false;
                    }
                    if xpp != NONE {
                        if let Some(l) = self.link_mut(xpp) {
                            l.red = true;
                        }
                        root = self.rotate_left(root, xpp);
                    }
                }
            }
        }
        root
    }

    /// Ports `TreeNode.removeTreeNode` (`HashMap.java:2190-2286`).
    fn remove_tree_node(&mut self, index: usize, p: u32, movable: bool) {
        if self.buckets.is_empty() {
            return;
        }
        let mut first = self.buckets.get(index).copied().unwrap_or(NONE);
        let mut root = first;
        let succ = self.nodes.get(p as usize).map_or(NONE, |nd| nd.next);
        let pred = self.link(p).map_or(NONE, |l| l.prev);

        if pred == NONE {
            first = succ;
            if let Some(b) = self.buckets.get_mut(index) {
                *b = first;
            }
        } else if let Some(pred_node) = self.nodes.get_mut(pred as usize) {
            pred_node.next = succ;
        }
        if succ != NONE {
            if let Some(succ_link) = self.link_mut(succ) {
                succ_link.prev = pred;
            }
        }
        if first == NONE {
            self.active_tree_bins = self.active_tree_bins.saturating_sub(1);
            self.reset_link(p);
            return;
        }
        if self.link(root).is_some_and(|l| l.parent != NONE) {
            root = self.tree_root(root);
        }
        let untreeify_now = root == NONE || {
            movable
                && self.link(root).is_none_or(|rl_root| {
                    rl_root.right == NONE
                        || rl_root.left == NONE
                        || self
                            .link(rl_root.left)
                            .is_none_or(|rl_left| rl_left.left == NONE)
                })
        };
        if untreeify_now {
            self.untreeify_chain(first);
            if let Some(b) = self.buckets.get_mut(index) {
                *b = first;
            }
            self.active_tree_bins = self.active_tree_bins.saturating_sub(1);
            self.reset_link(p);
            return;
        }

        let (pl, pr) = self.link(p).map_or((NONE, NONE), |l| (l.left, l.right));
        let replacement = if pl != NONE && pr != NONE {
            let mut s = pr;
            for _ in 0..=self.nodes.len() {
                let sl = self.link(s).map_or(NONE, |l| l.left);
                if sl == NONE {
                    break;
                }
                s = sl;
            }
            let s_red = self.link(s).is_some_and(|l| l.red);
            let p_red = self.link(p).is_some_and(|l| l.red);
            if let Some(slk) = self.link_mut(s) {
                slk.red = p_red;
            }
            if let Some(plk) = self.link_mut(p) {
                plk.red = s_red;
            }
            let sr = self.link(s).map_or(NONE, |l| l.right);
            let pp = self.link(p).map_or(NONE, |l| l.parent);
            if s == pr {
                if let Some(plk) = self.link_mut(p) {
                    plk.parent = s;
                }
                if let Some(slk) = self.link_mut(s) {
                    slk.right = p;
                }
            } else {
                let sp = self.link(s).map_or(NONE, |l| l.parent);
                if let Some(plk) = self.link_mut(p) {
                    plk.parent = sp;
                }
                if sp != NONE {
                    let s_was_left = self.link(sp).is_some_and(|l| s == l.left);
                    if let Some(splk) = self.link_mut(sp) {
                        if s_was_left {
                            splk.left = p;
                        } else {
                            splk.right = p;
                        }
                    }
                }
                if let Some(slk) = self.link_mut(s) {
                    slk.right = pr;
                }
                if let Some(prlk) = self.link_mut(pr) {
                    prlk.parent = s;
                }
            }
            if let Some(plk) = self.link_mut(p) {
                plk.left = NONE;
                plk.right = sr;
            }
            if sr != NONE {
                if let Some(srlk) = self.link_mut(sr) {
                    srlk.parent = p;
                }
            }
            if let Some(slk) = self.link_mut(s) {
                slk.left = pl;
                slk.parent = pp;
            }
            if let Some(pllk) = self.link_mut(pl) {
                pllk.parent = s;
            }
            if pp == NONE {
                root = s;
            } else {
                let p_was_left = self.link(pp).is_some_and(|l| p == l.left);
                if let Some(pplk) = self.link_mut(pp) {
                    if p_was_left {
                        pplk.left = s;
                    } else {
                        pplk.right = s;
                    }
                }
            }
            if sr != NONE {
                sr
            } else {
                p
            }
        } else if pl != NONE {
            pl
        } else if pr != NONE {
            pr
        } else {
            p
        };

        if replacement != p {
            let pp = self.link(p).map_or(NONE, |l| l.parent);
            if let Some(repl_lk) = self.link_mut(replacement) {
                repl_lk.parent = pp;
            }
            if pp == NONE {
                root = replacement;
                if let Some(repl_lk) = self.link_mut(replacement) {
                    repl_lk.red = false;
                }
            } else {
                let p_was_left = self.link(pp).is_some_and(|l| p == l.left);
                if let Some(pplk) = self.link_mut(pp) {
                    if p_was_left {
                        pplk.left = replacement;
                    } else {
                        pplk.right = replacement;
                    }
                }
            }
            if let Some(plk) = self.link_mut(p) {
                plk.left = NONE;
                plk.right = NONE;
                plk.parent = NONE;
            }
        }

        let p_is_red = self.link(p).is_some_and(|l| l.red);
        let r = if p_is_red {
            root
        } else {
            self.balance_deletion(root, replacement)
        };

        if replacement == p {
            let pp = self.link(p).map_or(NONE, |l| l.parent);
            if let Some(plk) = self.link_mut(p) {
                plk.parent = NONE;
            }
            if pp != NONE {
                let p_was_left = self.link(pp).is_some_and(|l| p == l.left);
                let p_was_right = self.link(pp).is_some_and(|l| p == l.right);
                if let Some(pplk) = self.link_mut(pp) {
                    if p_was_left {
                        pplk.left = NONE;
                    } else if p_was_right {
                        pplk.right = NONE;
                    }
                }
            }
        }
        if movable {
            self.move_root_to_front(index, r);
        }
        self.reset_link(p);
    }

    /// Ports `TreeNode.balanceDeletion` (`HashMap.java:2439-2529`).
    fn balance_deletion(&mut self, mut root: u32, mut x: u32) -> u32 {
        for _ in 0..=self.nodes.len().saturating_add(1) {
            if x == NONE || x == root {
                return root;
            }
            let xp = self.link(x).map_or(NONE, |l| l.parent);
            if xp == NONE {
                if let Some(xl) = self.link_mut(x) {
                    xl.red = false;
                }
                return x;
            }
            if self.link(x).is_some_and(|l| l.red) {
                if let Some(xl) = self.link_mut(x) {
                    xl.red = false;
                }
                return root;
            }
            let (mut xpl, mut xpr) = self.link(xp).map_or((NONE, NONE), |l| (l.left, l.right));
            if xpl == x {
                if xpr != NONE && self.link(xpr).is_some_and(|l| l.red) {
                    if let Some(l) = self.link_mut(xpr) {
                        l.red = false;
                    }
                    if let Some(l) = self.link_mut(xp) {
                        l.red = true;
                    }
                    root = self.rotate_left(root, xp);
                    let xp_now = self.link(x).map_or(NONE, |l| l.parent);
                    xpr = self.link(xp_now).map_or(NONE, |l| l.right);
                }
                let xp = self.link(x).map_or(NONE, |l| l.parent);
                if xpr == NONE {
                    x = xp;
                } else {
                    let (sl, sr) = self.link(xpr).map_or((NONE, NONE), |l| (l.left, l.right));
                    let sr_red = sr != NONE && self.link(sr).is_some_and(|l| l.red);
                    let sl_red = sl != NONE && self.link(sl).is_some_and(|l| l.red);
                    if !sr_red && !sl_red {
                        if let Some(l) = self.link_mut(xpr) {
                            l.red = true;
                        }
                        x = xp;
                    } else {
                        if !sr_red {
                            if sl != NONE {
                                if let Some(l) = self.link_mut(sl) {
                                    l.red = false;
                                }
                            }
                            if let Some(l) = self.link_mut(xpr) {
                                l.red = true;
                            }
                            root = self.rotate_right(root, xpr);
                            let xp_now = self.link(x).map_or(NONE, |l| l.parent);
                            xpr = self.link(xp_now).map_or(NONE, |l| l.right);
                        }
                        let xp = self.link(x).map_or(NONE, |l| l.parent);
                        if xpr != NONE {
                            let xp_red = xp != NONE && self.link(xp).is_some_and(|l| l.red);
                            if let Some(l) = self.link_mut(xpr) {
                                l.red = xp_red;
                            }
                            let sr_now = self.link(xpr).map_or(NONE, |l| l.right);
                            if sr_now != NONE {
                                if let Some(l) = self.link_mut(sr_now) {
                                    l.red = false;
                                }
                            }
                        }
                        if xp != NONE {
                            if let Some(l) = self.link_mut(xp) {
                                l.red = false;
                            }
                            root = self.rotate_left(root, xp);
                        }
                        x = root;
                    }
                }
            } else {
                if xpl != NONE && self.link(xpl).is_some_and(|l| l.red) {
                    if let Some(l) = self.link_mut(xpl) {
                        l.red = false;
                    }
                    if let Some(l) = self.link_mut(xp) {
                        l.red = true;
                    }
                    root = self.rotate_right(root, xp);
                    let xp_now = self.link(x).map_or(NONE, |l| l.parent);
                    xpl = self.link(xp_now).map_or(NONE, |l| l.left);
                }
                let xp = self.link(x).map_or(NONE, |l| l.parent);
                if xpl == NONE {
                    x = xp;
                } else {
                    let (sl, sr) = self.link(xpl).map_or((NONE, NONE), |l| (l.left, l.right));
                    let sl_red = sl != NONE && self.link(sl).is_some_and(|l| l.red);
                    let sr_red = sr != NONE && self.link(sr).is_some_and(|l| l.red);
                    if !sl_red && !sr_red {
                        if let Some(l) = self.link_mut(xpl) {
                            l.red = true;
                        }
                        x = xp;
                    } else {
                        if !sl_red {
                            if sr != NONE {
                                if let Some(l) = self.link_mut(sr) {
                                    l.red = false;
                                }
                            }
                            if let Some(l) = self.link_mut(xpl) {
                                l.red = true;
                            }
                            root = self.rotate_left(root, xpl);
                            let xp_now = self.link(x).map_or(NONE, |l| l.parent);
                            xpl = self.link(xp_now).map_or(NONE, |l| l.left);
                        }
                        let xp = self.link(x).map_or(NONE, |l| l.parent);
                        if xpl != NONE {
                            let xp_red = xp != NONE && self.link(xp).is_some_and(|l| l.red);
                            if let Some(l) = self.link_mut(xpl) {
                                l.red = xp_red;
                            }
                            let sl_now = self.link(xpl).map_or(NONE, |l| l.left);
                            if sl_now != NONE {
                                if let Some(l) = self.link_mut(sl_now) {
                                    l.red = false;
                                }
                            }
                        }
                        if xp != NONE {
                            if let Some(l) = self.link_mut(xp) {
                                l.red = false;
                            }
                            root = self.rotate_right(root, xp);
                        }
                        x = root;
                    }
                }
            }
        }
        root
    }

    /// Ports `TreeNode.split` (`HashMap.java:2298-2343`), updating `active_tree_bins`
    /// accurately across untreeifying, asymmetric, and symmetric splits.
    fn split_tree_bin(&mut self, index: usize, head: u32, old_cap: usize) {
        let mut lo_head = NONE;
        let mut lo_tail = NONE;
        let mut hi_head = NONE;
        let mut hi_tail = NONE;
        let mut lc = 0usize;
        let mut hc = 0usize;

        let mut e = head;
        for _ in 0..=self.nodes.len() {
            if e == NONE {
                break;
            }
            let (next, hash) = match self.nodes.get_mut(e as usize) {
                Some(nd) => {
                    let nxt = nd.next;
                    nd.next = NONE;
                    let h = Self::spread(&nd.key);
                    (nxt, h)
                }
                None => break,
            };
            if (hash & (old_cap as u32)) == 0 {
                if let Some(el) = self.link_mut(e) {
                    el.prev = lo_tail;
                }
                if lo_tail == NONE {
                    lo_head = e;
                } else if let Some(lt) = self.nodes.get_mut(lo_tail as usize) {
                    lt.next = e;
                }
                lo_tail = e;
                lc += 1;
            } else {
                if let Some(el) = self.link_mut(e) {
                    el.prev = hi_tail;
                }
                if hi_tail == NONE {
                    hi_head = e;
                } else if let Some(ht) = self.nodes.get_mut(hi_tail as usize) {
                    ht.next = e;
                }
                hi_tail = e;
                hc += 1;
            }
            e = next;
        }

        // Subtract the 1 pre-split tree bin at `index`; add back 1 for each half that stays a tree bin.
        self.active_tree_bins = self.active_tree_bins.saturating_sub(1);

        if lo_head != NONE {
            if lc <= UNTREEIFY_THRESHOLD {
                self.untreeify_chain(lo_head);
                if let Some(b) = self.buckets.get_mut(index) {
                    *b = lo_head;
                }
            } else {
                if let Some(b) = self.buckets.get_mut(index) {
                    *b = lo_head;
                }
                if hi_head != NONE {
                    self.treeify_bucket(index);
                }
                self.active_tree_bins = self.active_tree_bins.saturating_add(1);
            }
        }
        if hi_head != NONE {
            let hi_idx = index + old_cap;
            if hc <= UNTREEIFY_THRESHOLD {
                self.untreeify_chain(hi_head);
                if let Some(b) = self.buckets.get_mut(hi_idx) {
                    *b = hi_head;
                }
            } else {
                if let Some(b) = self.buckets.get_mut(hi_idx) {
                    *b = hi_head;
                }
                if lo_head != NONE {
                    self.treeify_bucket(hi_idx);
                }
                self.active_tree_bins = self.active_tree_bins.saturating_add(1);
            }
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

        let Some((new_cap, new_thr)) = next_resize_plan(old_cap) else {
            self.threshold = i32::MAX as usize;
            return;
        };
        self.threshold = new_thr;
        self.buckets.resize(new_cap, NONE);

        for j in 0..old_cap {
            let e = self.buckets.get(j).copied().unwrap_or(NONE);
            if e == NONE {
                continue;
            }
            if let Some(b) = self.buckets.get_mut(j) {
                *b = NONE;
            }
            let e_next = self.nodes.get(e as usize).map_or(NONE, |nd| nd.next);
            if e_next == NONE {
                let hash = self
                    .nodes
                    .get(e as usize)
                    .map_or(0, |nd| Self::spread(&nd.key));
                let dst = (hash as usize) & (new_cap - 1);
                if let Some(b) = self.buckets.get_mut(dst) {
                    *b = e;
                }
            } else if self.is_tree_node(e) {
                self.split_tree_bin(j, e, old_cap);
            } else {
                let mut curr = e;
                let mut lo_head = NONE;
                let mut lo_tail = NONE;
                let mut hi_head = NONE;
                let mut hi_tail = NONE;

                for _ in 0..=self.nodes.len() {
                    if curr == NONE {
                        break;
                    }
                    let Some(node) = self.nodes.get(curr as usize) else {
                        break;
                    };
                    let next = node.next;
                    let hash = Self::spread(&node.key);
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
                    if let Some(b) = self.buckets.get_mut(j) {
                        *b = lo_head;
                    }
                }
                if hi_tail != NONE {
                    if let Some(tail_node) = self.nodes.get_mut(hi_tail as usize) {
                        tail_node.next = NONE;
                    }
                    if let Some(b) = self.buckets.get_mut(j + old_cap) {
                        *b = hi_head;
                    }
                }
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

    pub(crate) fn into_entries(self) -> impl Iterator<Item = (K, V)> {
        IntoEntries {
            buckets: self.buckets,
            nodes: self.nodes,
            bucket_idx: 0,
            curr: NONE,
            remaining: self.size,
        }
    }
}

struct IntoEntries<K, V> {
    buckets: Vec<u32>,
    nodes: Vec<Node<K, V>>,
    bucket_idx: usize,
    curr: u32,
    remaining: usize,
}

impl<K, V> Iterator for IntoEntries<K, V>
where
    K: Clone,
    V: Default,
{
    type Item = (K, V);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        while self.curr == NONE {
            let &head = self.buckets.get(self.bucket_idx)?;
            self.bucket_idx += 1;
            self.curr = head;
        }
        let node = self.nodes.get_mut(self.curr as usize)?;
        self.curr = node.next;
        self.remaining -= 1;
        let k = node.key.clone();
        let v = std::mem::take(&mut node.val);
        Some((k, v))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<K, V> Default for JavaHashMap<K, V>
where
    K: Clone + JavaHash + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<K, V> IntoIterator for JavaHashMap<K, V>
where
    K: Clone + JavaHash + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    type Item = (K, V);
    type IntoIter = std::vec::IntoIter<(K, V)>;

    fn into_iter(mut self) -> Self::IntoIter {
        let mut items = Vec::with_capacity(self.size);
        let max_steps = self.nodes.len();
        for &head in &self.buckets {
            let mut curr = head;
            for _ in 0..=max_steps {
                if curr == NONE {
                    break;
                }
                let Some(node) = self.nodes.get_mut(curr as usize) else {
                    break;
                };
                let k = node.key.clone();
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
#[must_use]
pub fn java_string_hash(s: &str) -> i32 {
    if s.is_ascii() {
        return s
            .as_bytes()
            .iter()
            .fold(0i32, |h, &b| h.wrapping_mul(31).wrapping_add(i32::from(b)));
    }
    s.encode_utf16()
        .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(i32::from(c)))
}

/// Orders unique keys (given in first-insertion order) in OpenJDK 8+ `HashMap<String, V>`
/// `keySet()` / `entrySet()` iteration order, matching `HashMap.putVal` bucket chaining,
/// `treeifyBin` capacity-doubling when a bucket reaches 9 entries (`cap < 64`),
/// `TreeNode.treeify` / `putTreeVal` / `moveRootToFront` red-black tree ordering (`cap >= 64`),
/// and `loHead`/`hiHead` (`TreeNode.split`) on resize.
#[must_use]
pub fn java_hashmap_key_order<'a, S: AsRef<str> + 'a + ?Sized, I: IntoIterator<Item = &'a S>>(
    keys: I,
) -> Vec<&'a S> {
    let mut refs: Vec<&'a S> = Vec::new();
    let mut map: JavaHashMap<&'a str, usize> = JavaHashMap::new();
    for k in keys {
        let next_idx = refs.len();
        let stored = *map.get_or_put_with(k.as_ref(), || next_idx);
        if stored == next_idx {
            refs.push(k);
        }
    }
    map.iter()
        .filter_map(|(_, &i)| refs.get(i).copied())
        .collect()
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

        fn spread(key: &K) -> u32 {
            let h = (*key).into() as u32;
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
            let hash = Self::spread(key);
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
            let hash = Self::spread(key);
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
            let hash = Self::spread(&key);
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
            let hash = Self::spread(&key);
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

            let hash = Self::spread(&key);
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
            let hash = Self::spread(key);
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
                        let hash = Self::spread(&k);
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

    /// Verifies OpenJDK `TreeNode.checkInvariants` (`HashMap.java:2534-2554`), Red-Black
    /// black-height balance, `active_tree_bins` accuracy, and `free_head` `TreeLink` hygiene.
    fn assert_openjdk_invariants<K, V>(map: &JavaHashMap<K, V>)
    where
        K: Clone + JavaHash + Eq + PartialEq + Hash + Debug,
        V: Clone + Debug + Default,
    {
        if let Some(links) = map.tree_links.as_deref() {
            assert_eq!(
                links.len(),
                map.nodes.len(),
                "tree_links length must stay synchronized with nodes arena"
            );
        } else {
            assert_eq!(
                map.active_tree_bins, 0,
                "active_tree_bins must be 0 when tree_links is None"
            );
        }

        let mut free_count = 0usize;
        let mut fc = map.free_head;
        while fc != NONE {
            if let Some(links) = map.tree_links.as_deref() {
                assert_eq!(
                    links[fc as usize],
                    TreeLink::default(),
                    "free-list slot {fc} must have reset TreeLink::default()"
                );
            }
            fc = map.nodes[fc as usize].next;
            free_count += 1;
        }
        assert_eq!(map.len() + free_count, map.nodes.len());

        fn check_rb_subtree<K, V>(
            map: &JavaHashMap<K, V>,
            links: &[TreeLink],
            t: u32,
            visited: &mut std::collections::HashSet<u32>,
        ) -> usize
        where
            K: Clone + JavaHash + Eq + PartialEq + Hash + Debug,
            V: Clone + Debug + Default,
        {
            assert!(visited.insert(t), "cycle in TreeNode tree links at {t}");
            let lk = links[t as usize];
            assert!(lk.is_tree, "tree-reachable node {t} must have is_tree=true");
            let tn = map.nodes[t as usize].next;
            if lk.prev != NONE {
                assert_eq!(map.nodes[lk.prev as usize].next, t);
            }
            if tn != NONE {
                assert_eq!(links[tn as usize].prev, t);
            }
            if lk.parent != NONE {
                let plk = links[lk.parent as usize];
                assert!(plk.left == t || plk.right == t);
            }
            let th = JavaHashMap::<K, V>::spread(&map.nodes[t as usize].key) as i32;
            let bh_left = if lk.left != NONE {
                let llk = links[lk.left as usize];
                assert_eq!(llk.parent, t);
                let lh = JavaHashMap::<K, V>::spread(&map.nodes[lk.left as usize].key) as i32;
                assert!(lh <= th, "BST violation left hash {lh} > parent {th}");
                if lk.red {
                    assert!(!llk.red, "red node {t} has red left child {}", lk.left);
                }
                check_rb_subtree(map, links, lk.left, visited)
            } else {
                1
            };
            let bh_right = if lk.right != NONE {
                let rlk = links[lk.right as usize];
                assert_eq!(rlk.parent, t);
                let rh = JavaHashMap::<K, V>::spread(&map.nodes[lk.right as usize].key) as i32;
                assert!(rh >= th, "BST violation right hash {rh} < parent {th}");
                if lk.red {
                    assert!(!rlk.red, "red node {t} has red right child {}", lk.right);
                }
                check_rb_subtree(map, links, lk.right, visited)
            } else {
                1
            };
            assert_eq!(bh_left, bh_right, "Red-Black black-height mismatch at {t}");
            bh_left + usize::from(!lk.red)
        }

        let mut observed_tree_bins = 0u32;
        let mut total_live = 0usize;
        for &head in &map.buckets {
            if head == NONE {
                continue;
            }
            if map.is_tree_node(head) {
                observed_tree_bins += 1;
                let links = map
                    .tree_links
                    .as_deref()
                    .expect("tree bin requires tree_links");
                let root_lk = links[head as usize];
                assert_eq!(
                    root_lk.parent, NONE,
                    "bucket head must be tree root (moveRootToFront)"
                );
                assert_eq!(root_lk.prev, NONE, "bucket head must have prev == NONE");
                assert!(!root_lk.red, "Red-Black root must be black");

                let mut tree_set = std::collections::HashSet::new();
                let _ = check_rb_subtree(map, links, head, &mut tree_set);

                let mut list_set = std::collections::HashSet::new();
                let mut curr = head;
                while curr != NONE {
                    assert!(list_set.insert(curr));
                    curr = map.nodes[curr as usize].next;
                }
                assert_eq!(
                    tree_set, list_set,
                    "tree-reachable nodes must match next-chain nodes"
                );
                total_live += list_set.len();
            } else {
                let mut curr = head;
                while curr != NONE {
                    if let Some(links) = map.tree_links.as_deref() {
                        assert_eq!(
                            links[curr as usize],
                            TreeLink::default(),
                            "plain list node {curr} must have reset TreeLink::default()"
                        );
                    }
                    total_live += 1;
                    curr = map.nodes[curr as usize].next;
                }
            }
        }
        assert_eq!(observed_tree_bins, map.active_tree_bins);
        assert_eq!(total_live, map.len());
    }

    #[derive(Debug, PartialEq, Eq)]
    struct LiveJavaDump {
        len: usize,
        cap: usize,
        has_tree: bool,
        entries: Vec<(i32, i32)>,
    }

    fn find_java_toolchain() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let check_dir =
            |dir: std::path::PathBuf| -> Option<(std::path::PathBuf, std::path::PathBuf)> {
                let javac = dir.join("javac");
                let java = dir.join("java");
                if javac.is_file() && java.is_file() {
                    Some((javac, java))
                } else {
                    None
                }
            };

        if let Ok(java_home) = std::env::var("JAVA_HOME") {
            if let Some(pair) = check_dir(std::path::PathBuf::from(java_home).join("bin")) {
                return Some(pair);
            }
        }
        if let Ok(home) = std::env::var("HOME") {
            if let Some(pair) =
                check_dir(std::path::PathBuf::from(home).join("miniconda3/envs/vardict-oracle/bin"))
            {
                return Some(pair);
            }
        }
        if let Ok(path_var) = std::env::var("PATH") {
            for dir in std::env::split_paths(&path_var) {
                if let Some(pair) = check_dir(dir) {
                    return Some(pair);
                }
            }
        }
        None
    }

    fn run_live_java_oracle(batch_input: &str) -> Option<Vec<LiveJavaDump>> {
        use std::io::Write;
        use std::process::{Command, Stdio};

        let (javac, java) = find_java_toolchain()?;

        struct TempDirGuard(std::path::PathBuf);
        impl Drop for TempDirGuard {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }

        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "vardictrust_javamap_oracle_{}_{}",
            std::process::id(),
            id
        ));
        std::fs::create_dir_all(&dir).ok()?;
        let _guard = TempDirGuard(dir.clone());

        let java_src = r#"import java.io.*;
import java.lang.reflect.Field;
import java.util.*;

public class BatchJavaMapOracle {
    private static final Field TABLE_FIELD;
    static {
        try {
            TABLE_FIELD = HashMap.class.getDeclaredField("table");
            TABLE_FIELD.setAccessible(true);
        } catch (ReflectiveOperationException e) {
            throw new ExceptionInInitializerError(e);
        }
    }
    private static int cap(HashMap<Integer, Integer> m) {
        try {
            Object[] t = (Object[]) TABLE_FIELD.get(m);
            return t == null ? 0 : t.length;
        } catch (IllegalAccessException e) { throw new RuntimeException(e); }
    }
    private static boolean hasTree(HashMap<Integer, Integer> m) {
        try {
            Object[] t = (Object[]) TABLE_FIELD.get(m);
            if (t == null) return false;
            for (Object n : t) {
                if (n != null && n.getClass().getSimpleName().contains("TreeNode")) return true;
            }
            return false;
        } catch (IllegalAccessException e) { throw new RuntimeException(e); }
    }
    public static void main(String[] args) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        StringBuilder out = new StringBuilder();
        HashMap<Integer, Integer> m = new HashMap<>();
        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) continue;
            char c = line.charAt(0);
            if (c == 'N') {
                m = new HashMap<>();
            } else if (c == 'P') {
                int sp = line.indexOf(' ', 2);
                int k = Integer.parseInt(line.substring(2, sp));
                int v = Integer.parseInt(line.substring(sp + 1));
                m.put(k, v);
            } else if (c == 'G') {
                int sp1 = line.indexOf(' ', 2);
                int sp2 = line.indexOf(' ', sp1 + 1);
                int k = Integer.parseInt(line.substring(2, sp1));
                int v = Integer.parseInt(line.substring(sp1 + 1, sp2));
                int delta = Integer.parseInt(line.substring(sp2 + 1));
                Integer cur = m.get(k);
                if (cur == null) { cur = v; m.put(k, cur); }
                m.put(k, cur + delta);
            } else if (c == 'C') {
                int sp1 = line.indexOf(' ', 2);
                int sp2 = line.indexOf(' ', sp1 + 1);
                int k = Integer.parseInt(line.substring(2, sp1));
                int init = Integer.parseInt(line.substring(sp1 + 1, sp2));
                int delta = Integer.parseInt(line.substring(sp2 + 1));
                int cur = m.computeIfAbsent(k, x -> x + init);
                m.put(k, cur + delta);
            } else if (c == 'M') {
                int sp = line.indexOf(' ', 2);
                int k = Integer.parseInt(line.substring(2, sp));
                int delta = Integer.parseInt(line.substring(sp + 1));
                Integer cur = m.get(k);
                if (cur != null) m.put(k, cur + delta);
            } else if (c == 'R') {
                int k = Integer.parseInt(line.substring(2));
                m.remove(k);
            } else if (c == 'D') {
                out.append(m.size()).append(' ')
                   .append(cap(m)).append(' ')
                   .append(hasTree(m) ? '1' : '0');
                for (Map.Entry<Integer, Integer> e : m.entrySet()) {
                    out.append(' ').append(e.getKey()).append(':').append(e.getValue());
                }
                out.append('\n');
            }
        }
        System.out.print(out);
    }
}
"#;
        let src_path = dir.join("BatchJavaMapOracle.java");
        std::fs::write(&src_path, java_src).expect("write BatchJavaMapOracle.java");
        let compile = Command::new(&javac)
            .arg(&src_path)
            .output()
            .expect("spawn javac");
        assert!(
            compile.status.success(),
            "javac BatchJavaMapOracle failed: {}",
            String::from_utf8_lossy(&compile.stderr)
        );

        let ver_out = Command::new(&java)
            .arg("-version")
            .output()
            .expect("java -version");
        let ver_str = format!(
            "{}{}",
            String::from_utf8_lossy(&ver_out.stdout),
            String::from_utf8_lossy(&ver_out.stderr)
        );
        let mut cmd = Command::new(&java);
        if !ver_str.contains("version \"1.8.") {
            cmd.args(["--add-opens", "java.base/java.util=ALL-UNNAMED"]);
        }
        let mut child = cmd
            .args(["-cp", &dir.to_string_lossy(), "BatchJavaMapOracle"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("spawn java BatchJavaMapOracle");
        if let Some(mut stdin) = child.stdin.take() {
            stdin
                .write_all(batch_input.as_bytes())
                .expect("write BatchJavaMapOracle stdin");
        }
        let output = child
            .wait_with_output()
            .expect("wait BatchJavaMapOracle output");
        assert!(
            output.status.success(),
            "java BatchJavaMapOracle failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8(output.stdout).expect("utf8 BatchJavaMapOracle stdout");
        let mut dumps = Vec::new();
        for line in stdout.lines() {
            let mut parts = line.split_whitespace();
            let len: usize = parts.next().expect("len").parse().expect("parse len");
            let cap: usize = parts.next().expect("cap").parse().expect("parse cap");
            let has_tree = parts.next().expect("has_tree") == "1";
            let mut entries = Vec::with_capacity(len);
            for kv in parts {
                let (ks, vs) = kv.split_once(':').expect("k:v");
                entries.push((ks.parse().expect("parse k"), vs.parse().expect("parse v")));
            }
            dumps.push(LiveJavaDump {
                len,
                cap,
                has_tree,
                entries,
            });
        }
        Some(dumps)
    }

    #[test]
    #[expect(
        deprecated,
        reason = "verifies OpenJDK HashMap.computeIfAbsent C-op oracle streams and 1.0.x SemVer shim"
    )]
    fn javamap_arena_matches_legacy_on_random_op_streams() {
        use std::fmt::Write as _;
        fn next_u64(state: &mut u64) -> u64 {
            let mut x = *state;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            *state = x;
            x
        }

        // First 20 streams stay strictly below treeify threshold -> compare step-by-step against JavaHashMapLegacy!
        for stream_idx in 0..20u64 {
            let mut rng =
                0x517C_C1B7_2722_0A95u64 ^ (stream_idx.wrapping_mul(0xBF58_476D_1CE4_E5B9));
            if rng == 0 {
                rng = 1;
            }
            let mut arena: JavaHashMap<i32, i32> = JavaHashMap::new();
            let mut legacy: JavaHashMapLegacy<i32, i32> = JavaHashMapLegacy::new();
            for step in 0..2_000usize {
                let r = next_u64(&mut rng);
                let op = r % 100;
                // Consecutive small integer keys distribute across buckets (<= 2 per bucket, never treeifying)
                let key = ((r >> 12) % 96) as i32 - 32;
                let val = (step as i32).wrapping_mul(3) ^ key;
                if op < 30 {
                    assert_eq!(arena.insert(key, val), legacy.insert(key, val));
                } else if op < 60 {
                    let delta = ((r >> 32) & 0xFF) as i32;
                    let a_ref = arena.compute_if_absent(key, |k| k.wrapping_add(val));
                    *a_ref = a_ref.wrapping_add(delta);
                    let l_ref = legacy.compute_if_absent(key, |k| k.wrapping_add(val));
                    *l_ref = l_ref.wrapping_add(delta);
                    assert_eq!(*a_ref, *l_ref);
                } else if op < 80 {
                    assert_eq!(arena.remove(&key), legacy.remove(&key));
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
                assert!(!arena.treeified() && !legacy.treeified());
                assert!(
                    arena.tree_links.is_none(),
                    "non-treeifying maps must never allocate tree_links"
                );
                if step % 500 == 0 || step == 1_999 {
                    assert_openjdk_invariants(&arena);
                    let a_iter: Vec<(i32, i32)> = arena.iter().map(|(&k, &v)| (k, v)).collect();
                    let l_iter: Vec<(i32, i32)> = legacy.iter().map(|(&k, &v)| (k, v)).collect();
                    assert_eq!(a_iter, l_iter);
                }
            }
        }

        // Next 50 heavy-collision streams exercise TreeNode treeify/split/untreeify ->
        // compare non-treeified prefix against JavaHashMapLegacy and all checkpoints against OpenJDK JavaMapOracle!
        let mut oracle_script = String::new();
        let mut arena_snapshots: Vec<(u64, usize, LiveJavaDump)> = Vec::new();

        for stream_idx in 0..50u64 {
            let mut rng =
                0x9E37_79B9_7F4A_7C15u64 ^ (stream_idx.wrapping_mul(0xBF58_476D_1CE4_E5B9));
            if rng == 0 {
                rng = 1;
            }
            let mut arena: JavaHashMap<i32, i32> = JavaHashMap::new();
            let mut legacy: JavaHashMapLegacy<i32, i32> = JavaHashMapLegacy::new();
            oracle_script.push_str("N\n");

            for step in 0..5_000usize {
                let r = next_u64(&mut rng);
                let op = r % 100;
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
                    let _ = writeln!(oracle_script, "P {key} {val}");
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
                    let _ = writeln!(oracle_script, "C {key} {val} {delta}");
                    assert_eq!(
                        *a_ref, *l_ref,
                        "compute_if_absent value mismatch at stream {stream_idx} step {step}"
                    );
                } else if op < 80 {
                    let a_rem = arena.remove(&key);
                    let l_rem = legacy.remove(&key);
                    let _ = writeln!(oracle_script, "R {key}");
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
                    let _ = writeln!(oracle_script, "M {key} 7");
                } else {
                    arena = arena.clone();
                    legacy = legacy.clone();
                }

                if step % 500 == 0 || step == 4_999 {
                    assert_openjdk_invariants(&arena);
                    assert_eq!(arena.len(), legacy.len());
                    assert_eq!(arena.is_empty(), legacy.is_empty());
                    let a_iter: Vec<(i32, i32)> = arena.iter().map(|(&k, &v)| (k, v)).collect();
                    let a_into: Vec<(i32, i32)> = arena.clone().into_iter().collect();
                    assert_eq!(a_iter, a_into);

                    if !legacy.treeified() {
                        assert!(
                            !arena.treeified(),
                            "non-treeifying prefix must keep treeified()==false at stream {stream_idx} step {step}"
                        );
                        let l_iter: Vec<(i32, i32)> =
                            legacy.iter().map(|(&k, &v)| (k, v)).collect();
                        assert_eq!(
                            a_iter, l_iter,
                            "non-treeifying iter() mismatch vs JavaHashMapLegacy at stream {stream_idx} step {step}"
                        );
                    }

                    oracle_script.push_str("D\n");
                    arena_snapshots.push((
                        stream_idx,
                        step,
                        LiveJavaDump {
                            len: arena.len(),
                            cap: arena.capacity(),
                            has_tree: arena.treeified(),
                            entries: a_iter,
                        },
                    ));
                }
            }
        }

        if let Some(java_dumps) = run_live_java_oracle(&oracle_script) {
            assert_eq!(java_dumps.len(), arena_snapshots.len());
            for ((stream_idx, step, arena_dump), java_dump) in
                arena_snapshots.iter().zip(java_dumps.iter())
            {
                assert_eq!(
                    arena_dump, java_dump,
                    "OpenJDK JavaMapOracle mismatch at stream {stream_idx} step {step}"
                );
            }
        }
    }

    #[test]
    #[expect(
        deprecated,
        reason = "replays OpenJDK HashMap.computeIfAbsent C-op traces in ops_1/ops_2/ops_treeify"
    )]
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
        let order_treeify = ops_treeify
            .lines()
            .last()
            .and_then(|l| l.strip_prefix("ORDER "));
        assert!(order_treeify.is_some(), "ops_treeify must end with ORDER");
        let m = run_oracle(ops_treeify, order_treeify);
        assert!(m.treeified(), "ops_treeify must set treeified flag");
    }

    #[expect(
        deprecated,
        reason = "replays OpenJDK HashMap.computeIfAbsent C-op traces in ops_put_threshold/ops_getput"
    )]
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
                }
                assert_eq!(m.capacity(), expected_cap, "capacity mismatch on {label}");
                let actual_order = m
                    .keys()
                    .map(|k| k.to_string())
                    .collect::<Vec<_>>()
                    .join(" ");
                assert_eq!(
                    actual_order, expected_order,
                    "key order mismatch on {label} (regime={regime})"
                );
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
    #[expect(
        deprecated,
        reason = "verifies OpenJDK HashMap.computeIfAbsent 8-node treeify threshold"
    )]
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
        use std::fmt::Write as _;
        fn next_u64(state: &mut u64) -> u64 {
            let mut x = *state;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            *state = x;
            x
        }

        // Non-treeifying streams compared directly against JavaHashMapLegacy:
        for stream_idx in 0..20u64 {
            let mut rng =
                0xC3C3_3C3C_5A5A_A5A5u64 ^ (stream_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15));
            if rng == 0 {
                rng = 1;
            }
            let mut arena: JavaHashMap<i32, i32> = JavaHashMap::new();
            let mut legacy: JavaHashMapLegacy<i32, i32> = JavaHashMapLegacy::new();
            for step in 0..1_000usize {
                let r = next_u64(&mut rng);
                let key = ((r >> 10) % 120) as i32;
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
                assert!(!arena.treeified() && !legacy.treeified());
                if step % 200 == 0 || step == 999 {
                    let a_into: Vec<(i32, i32)> = arena.clone().into_iter().collect();
                    let l_into: Vec<(i32, i32)> = legacy.clone().into_vec();
                    assert_eq!(a_into, l_into);
                }
            }
        }

        // Treeifying streams: non-treeifying prefix matches JavaHashMapLegacy; all checkpoints match OpenJDK JavaMapOracle!
        let mut oracle_script = String::new();
        let mut arena_snapshots: Vec<(u64, usize, LiveJavaDump)> = Vec::new();

        for stream_idx in 0..40u64 {
            let mut rng =
                0xA5A5_5A5A_F0F0_0F0Fu64 ^ (stream_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15));
            if rng == 0 {
                rng = 1;
            }
            let mut arena: JavaHashMap<i32, i32> = JavaHashMap::new();
            let mut legacy: JavaHashMapLegacy<i32, i32> = JavaHashMapLegacy::new();
            oracle_script.push_str("N\n");

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
                    let _ = writeln!(oracle_script, "G {key} {val} {delta}");
                    assert_eq!(*a_ref, *l_ref);
                } else {
                    assert_eq!(arena.remove(&key), legacy.remove(&key));
                    let _ = writeln!(oracle_script, "R {key}");
                }
                assert_eq!(arena.len(), legacy.len());
                if !legacy.treeified() {
                    assert!(!arena.treeified());
                }
                if step % 200 == 0 || step == 1_999 {
                    assert_openjdk_invariants(&arena);
                    let a_into: Vec<(i32, i32)> = arena.clone().into_iter().collect();
                    if !legacy.treeified() {
                        let l_into: Vec<(i32, i32)> = legacy.clone().into_vec();
                        assert_eq!(a_into, l_into);
                    }
                    oracle_script.push_str("D\n");
                    arena_snapshots.push((
                        stream_idx,
                        step,
                        LiveJavaDump {
                            len: arena.len(),
                            cap: arena.capacity(),
                            has_tree: arena.treeified(),
                            entries: a_into,
                        },
                    ));
                }
            }
        }

        if let Some(java_dumps) = run_live_java_oracle(&oracle_script) {
            assert_eq!(java_dumps.len(), arena_snapshots.len());
            for ((stream_idx, step, arena_dump), java_dump) in
                arena_snapshots.iter().zip(java_dumps.iter())
            {
                assert_eq!(
                    arena_dump, java_dump,
                    "OpenJDK JavaMapOracle get_or_put mismatch at stream {stream_idx} step {step}"
                );
            }
        }
    }

    #[test]
    #[expect(
        deprecated,
        reason = "replays OpenJDK HashMap.computeIfAbsent C-op trace in ops_treeify"
    )]
    fn op_log_matches_java_oracle_treeify() {
        let ops_treeify = include_str!("../testdata/javamap/ops_treeify.txt");
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        let mut expected_order = None;
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
            } else if let Some(ord) = line.strip_prefix("ORDER ") {
                expected_order = Some(ord);
            }
        }
        assert!(m.treeified(), "ops_treeify must set treeified flag");
        let actual = m
            .keys()
            .map(|k| k.to_string())
            .collect::<Vec<_>>()
            .join(" ");
        assert_eq!(
            Some(actual.as_str()),
            expected_order,
            "ops_treeify ORDER mismatch"
        );
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

        // Opt-in new_key_log captures absent-key allocations across insert / get_or_put_with
        // (including free_head reuse and across resize), ignores present-key updates, and does
        // not propagate across Clone.
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
        logged.get_or_put_with(30, || 88); // existing -> not logged
        logged.get_or_put_with(40, || 5); // new -> logged
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
        map.get_or_put_with(8 * 1024 + 7, || 8 * 1024 + 7);
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

    #[test]
    fn batch12_maximum_capacity_clamp_matches_openjdk() {
        // Simulate a map already at OpenJDK MAXIMUM_CAPACITY = 1 << 30 without allocating 4 GiB:
        // calling resize when old_cap >= 1 << 30 must clamp threshold to i32::MAX and return
        // without doubling old_cap << 1.
        let mut map: JavaHashMap<i32, i32> = JavaHashMap::new();
        map.insert(1, 10);
        assert_eq!(map.capacity(), 16);
        assert_eq!(map.threshold, 12);
        let cap_at_max = 1usize << 30;
        let expected_openjdk_thr = i32::MAX as usize;
        assert_eq!(
            super::next_resize_plan(1usize << 29),
            Some((cap_at_max, expected_openjdk_thr))
        );
        assert_eq!(super::next_resize_plan(cap_at_max), None);
    }

    #[test]
    fn tree_bin_signed_i32_hash_order_matches_openjdk() {
        // Bucket 63 at cap 64 with -64 and positive/negative colliding keys:
        let mut m63: JavaHashMap<i32, i32> = JavaHashMap::new();
        let k63 = [63, 127, 191, 255, 319, 383, -64, -128, -192, -256, -320];
        for (idx, &k) in k63.iter().enumerate() {
            assert_eq!(m63.insert(k, idx as i32), None);
        }
        assert!(m63.treeified());
        for (idx, &k) in k63.iter().enumerate() {
            assert_eq!(m63.get(&k), Some(&(idx as i32)));
        }
        assert_eq!(
            m63.keys().copied().collect::<Vec<_>>(),
            vec![127, 63, 191, 255, 319, 383, -64, -128, -192, -256, -320]
        );

        // Bucket 0 at cap 64 with -1 (spread == -65536) and non-negative multiples of 64:
        let mut m0: JavaHashMap<i32, i32> = JavaHashMap::new();
        let k0 = [0, 64, 128, 192, 256, 320, 384, 448, 512, 576, -1];
        for (idx, &k) in k0.iter().enumerate() {
            assert_eq!(m0.insert(k, idx as i32), None);
        }
        assert!(m0.treeified());
        assert_eq!(m0.get(&-1), Some(&10));
        assert_eq!(
            m0.keys().copied().collect::<Vec<_>>(),
            vec![192, 0, 64, 128, 256, 320, 384, 448, 512, 576, -1]
        );

        // Raw `(i << 16) | 0x8000_1234` bit-31 keys mixed with `(i << 16) | 0x0000_1234`:
        let mut m_raw: JavaHashMap<i32, i32> = JavaHashMap::new();
        for j in 0..6i32 {
            let pos = ((j * 64) << 16) | 0x0000_1234;
            m_raw.insert(pos, j);
        }
        for j in 0..6i32 {
            let neg = (((j * 64) as u32) << 16) | 0x8000_1234u32;
            m_raw.insert(neg as i32, 6 + j);
        }
        assert!(m_raw.treeified());
        assert_eq!(
            m_raw.keys().copied().collect::<Vec<_>>(),
            vec![
                4_198_964,
                4_660,
                8_393_268,
                12_587_572,
                16_781_876,
                20_976_180,
                -2_147_478_988,
                -2_143_284_684,
                -2_139_090_380,
                -2_134_896_076,
                -2_130_701_772,
                -2_126_507_468,
            ]
        );
    }

    #[test]
    #[expect(
        deprecated,
        reason = "verifies OpenJDK HashMap.computeIfAbsent tree-bin pre-resize and putTreeVal parity"
    )]
    fn tree_bin_compute_if_absent_matches_openjdk() {
        let mut cia: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..10 {
            let key = i * 64;
            let v = *cia.compute_if_absent(key, |&k| k + 1);
            assert_eq!(v, key + 1);
        }
        assert_eq!(cia.capacity(), 64);
        assert!(cia.treeified());
        assert_eq!(
            cia.keys().copied().collect::<Vec<_>>(),
            vec![384, 576, 512, 448, 320, 256, 192, 128, 64, 0]
        );

        // Grow to cap=64 with 39 distinct keys (1..=39), then treeify bin 0 on 8th compute_if_absent,
        // insert 9th & 10th via TreeNode.putTreeVal up to size == 49 > threshold (48) WITHOUT post-insert resize,
        // and verify the NEXT compute_if_absent pre-resizes 64 -> 128 and untreeifies.
        let mut cia64: JavaHashMap<i32, i32> = JavaHashMap::new();
        for k in 1..=39 {
            cia64.insert(k, k);
        }
        assert_eq!(cia64.capacity(), 64);
        assert!(!cia64.treeified());
        for i in 0..7 {
            let key = (10 - i) * 64;
            cia64.compute_if_absent(key, |&k| k * 2);
        }
        assert!(!cia64.treeified(), "7 colliding nodes must not treeify");
        cia64.compute_if_absent(3 * 64, |&k| k * 2); // 8th node in bin 0 -> treeifyBin!
        assert!(
            cia64.treeified(),
            "8th compute_if_absent node at cap 64 must treeify"
        );

        // Present lookup in tree bin must NOT invoke closure:
        let mut called = false;
        let existing = *cia64.compute_if_absent(3 * 64, |_| {
            called = true;
            -1
        });
        assert!(!called);
        assert_eq!(existing, 384);

        // Absent insertions into already-treeified bin 0 via putTreeVal:
        cia64.compute_if_absent(11 * 64, |&k| k * 2); // size == 48
        cia64.compute_if_absent(64, |&k| k * 2); // size == 49 > 48, no post-insert resize!
        assert_eq!(cia64.len(), 49);
        assert_eq!(cia64.capacity(), 64);
        assert!(cia64.treeified());
        let bin0_keys: Vec<i32> = cia64.keys().copied().filter(|&k| (k & 63) == 0).collect();
        assert_eq!(
            bin0_keys,
            vec![384, 192, 64, 256, 320, 448, 512, 576, 640, 704]
        );

        // Next compute_if_absent (even on an existing key) triggers pre-resize 64 -> 128,
        // splitting 10 keys into 5 lo + 5 hi (both <= UNTREEIFY_THRESHOLD=6 -> untreeified!).
        let hit = *cia64.compute_if_absent(64, |_| -999);
        assert_eq!(hit, 128);
        assert_eq!(cia64.capacity(), 128);
        assert!(!cia64.treeified());
    }

    #[test]
    fn tree_bin_cap64_49th_insert_treeifies_then_splits_matches_openjdk() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for k in 1..=40 {
            m.insert(k, k);
        }
        let nine = [64, 192, 320, 448, 576, 128, 256, 384, 512];
        for &k in &nine[..8] {
            m.insert(k, k * 10);
        }
        assert_eq!(m.len(), 48);
        assert_eq!(m.capacity(), 64);
        assert!(!m.treeified());

        // 49th insert is 9th node in bin 0: treeifyBin runs FIRST (moving root 192 ahead of 64),
        // then ++size > 48 immediately calls resize() 64 -> 128, splitting into 4 lo + 5 hi (untreeified).
        m.insert(nine[8], nine[8] * 10);
        assert_eq!(m.len(), 49);
        assert_eq!(m.capacity(), 128);
        assert!(!m.treeified());

        let b0: Vec<i32> = m.keys().copied().filter(|&k| (k & 127) == 0).collect();
        let b64: Vec<i32> = m.keys().copied().filter(|&k| (k & 127) == 64).collect();
        assert_eq!(b0, vec![128, 256, 384, 512]);
        assert_eq!(b64, vec![192, 64, 320, 448, 576]);
    }

    #[test]
    fn remove_tree_node_topology_untreeify_and_asymmetric_split_matches_openjdk() {
        // 1. Pointer-topology untreeification + 2-child interior root swap:
        let mut rm: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..8 {
            rm.insert(i * 128, i);
        }
        for i in 0..3 {
            rm.insert(i * 128 + 64, 100 + i);
        }
        assert_eq!(rm.len(), 11);
        assert_eq!(rm.capacity(), 64);
        assert!(rm.treeified());
        assert_eq!(
            rm.keys().copied().collect::<Vec<_>>(),
            vec![384, 0, 128, 256, 512, 640, 768, 896, 64, 192, 320]
        );

        // Remove root 384 (2-child interior node -> swaps tree links with successor 512 while preserving prev/next):
        assert_eq!(rm.remove(&384), Some(3));
        assert!(rm.treeified());
        assert_eq!(
            rm.keys().copied().collect::<Vec<_>>(),
            vec![512, 0, 128, 256, 640, 768, 896, 64, 192, 320]
        );

        assert_eq!(rm.remove(&256), Some(2));
        assert_eq!(rm.remove(&512), Some(4));
        assert_eq!(rm.remove(&768), Some(6));
        assert_eq!(rm.remove(&128), Some(1));
        // At 6 nodes, pointer topology STILL has root.right != NONE && root.left != NONE && root.left.left != NONE!
        assert_eq!(rm.len(), 6);
        assert!(
            rm.treeified(),
            "6-node tree with deep left grandchild must remain treeified per OpenJDK removeTreeNode"
        );
        assert_eq!(
            rm.keys().copied().collect::<Vec<_>>(),
            vec![640, 0, 896, 64, 192, 320]
        );

        // Even at 5 nodes after removing 0, topology still satisfies untreeify guard -> stays treeified!
        assert_eq!(rm.remove(&0), Some(0));
        assert_eq!(rm.len(), 5);
        assert!(
            rm.treeified(),
            "5-node tree with root.left.left != NONE must remain treeified per OpenJDK removeTreeNode"
        );
        assert_eq!(
            rm.keys().copied().collect::<Vec<_>>(),
            vec![640, 896, 64, 192, 320]
        );

        // Removing 64: BEFORE removing 64, root.left.left WAS 64 (!= NONE), so OpenJDK
        // checks topology prior to unlink and keeps the 4-node bin treeified!
        assert_eq!(rm.remove(&64), Some(100));
        assert_eq!(rm.len(), 4);
        assert!(
            rm.treeified(),
            "4-node bin stays treeified when root.left.left was non-null on removeTreeNode entry"
        );
        assert_eq!(
            rm.keys().copied().collect::<Vec<_>>(),
            vec![640, 896, 192, 320]
        );
        assert_openjdk_invariants(&rm);

        // Next removal (320): on entry root.left.left IS NONE -> OpenJDK untreeifies to plain list!
        assert_eq!(rm.remove(&320), Some(102));
        assert_eq!(rm.len(), 3);
        assert!(
            !rm.treeified(),
            "removing 320 triggers pointer-topology untreeify"
        );
        assert_eq!(rm.keys().copied().collect::<Vec<_>>(), vec![640, 896, 192]);
        assert_openjdk_invariants(&rm);

        // Free-list reuse after untreeify must allocate cleanly without stale tree flags:
        rm.insert(1, 777);
        assert!(!rm.treeified());
        assert_eq!(rm.get(&1), Some(&777));
        assert_openjdk_invariants(&rm);

        // 2. Asymmetric resize split (lc = 8 > 6 stays tree, hc = 3 <= 6 untreeifies -> active_tree_bins == 1!):
        let mut asym: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..8 {
            asym.insert(i * 128, i);
        }
        for i in 0..3 {
            asym.insert(i * 128 + 64, 100 + i);
        }
        assert!(asym.treeified());
        for k in 1..=38 {
            asym.insert(k, k);
        }
        assert_eq!(asym.len(), 49);
        assert_eq!(asym.capacity(), 128);
        assert!(
            asym.treeified(),
            "asymmetric split (lc=8 > 6, hc=3 <= 6) must retain active_tree_bins == 1"
        );
        assert_openjdk_invariants(&asym);
        let asym_b0: Vec<i32> = asym.keys().copied().filter(|&k| (k & 127) == 0).collect();
        let asym_b64: Vec<i32> = asym.keys().copied().filter(|&k| (k & 127) == 64).collect();
        assert_eq!(asym_b0, vec![384, 0, 128, 256, 512, 640, 768, 896]);
        assert_eq!(asym_b64, vec![64, 192, 320]);

        // 3. Symmetric resize split (lc = 8 > 6 and hc = 8 > 6 -> active_tree_bins == 2!):
        let mut sym: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..8 {
            sym.insert(i * 128, i);
            sym.insert(i * 128 + 64, 100 + i);
        }
        assert_eq!(sym.active_tree_bins, 1);
        for k in 1..=33 {
            sym.insert(k, k);
        }
        assert_eq!(sym.len(), 49);
        assert_eq!(sym.capacity(), 128);
        assert_eq!(sym.active_tree_bins, 2);
        assert_openjdk_invariants(&sym);

        // 4. Compact Node<i32, i32> size (zero extra bytes per node vs non-tree arena):
        assert_eq!(std::mem::size_of::<Node<i32, i32>>(), 12);
    }

    #[test]
    fn contains_both_and_first_precedes_after_move_root_to_front() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        // Insert 11 keys colliding in bin 0 at cap 64; 0 is inserted first, 384 becomes RB root & moves to front!
        for i in 0..11 {
            m.insert(i * 64, i);
        }
        assert!(m.treeified());
        let order: Vec<i32> = m.keys().copied().collect();
        assert_eq!(order.first().copied(), Some(192));
        assert!(m.contains_both_and_first_precedes(&192, &0));
        assert!(!m.contains_both_and_first_precedes(&0, &192));

        let rank: std::collections::HashMap<i32, usize> = order
            .iter()
            .copied()
            .enumerate()
            .map(|(i, k)| (k, i))
            .collect();
        let probes: Vec<i32> = (0..=13).map(|i| i * 64).collect();
        for &a in &probes {
            for &b in &probes {
                let expected =
                    matches!((rank.get(&a), rank.get(&b)), (Some(ia), Some(ib)) if ia < ib);
                assert_eq!(
                    m.contains_both_and_first_precedes(&a, &b),
                    expected,
                    "tree bin precedes mismatch for ({a}, {b})"
                );
            }
        }
    }

    #[test]
    fn tree_bin_utf16_non_bmp_string_order_matches_openjdk() {
        // Direct UTF-16 vs UTF-8 code-unit inversion check:
        // BMP U+E000 encodes to [0xE000]; non-BMP U+10000 encodes to surrogate pair [0xD800, 0xDC00].
        // In UTF-8 '\u{E000}' < '\u{10000}', whereas in Java String.compareTo (UTF-16) 0xE000 > 0xD800!
        assert_eq!(
            "\u{E000}".java_compare("\u{10000}"),
            std::cmp::Ordering::Greater
        );
        assert_eq!(
            "\u{10000}".java_compare("\u{E000}"),
            std::cmp::Ordering::Less
        );

        // Adversarial exact-hash collision in a TreeNode bin:
        // "\u{E000}Dk1PZv5" and "\u{10000}tUJg1oO" both have java_string_hash == -841_374_349!
        // Appending 6 identical-hash "Aa"/"BB" suffixes yields 12 distinct keys sharing the EXACT
        // same 32-bit hash, forcing TreeNode.treeify / putTreeVal / getTreeNode to call java_compare!
        let bmp_base = "\u{E000}Dk1PZv5";
        let supp_base = "\u{10000}tUJg1oO";
        assert_eq!(java_string_hash(bmp_base), -841_374_349);
        assert_eq!(java_string_hash(supp_base), -841_374_349);

        let suffixes = ["AaAaAa", "AaAaBB", "AaBBAa", "AaBBBB", "BBAaAa", "BBAaBB"];
        let mut m: JavaHashMap<String, i32> = JavaHashMap::new();
        for (idx, sfx) in suffixes.iter().enumerate() {
            m.insert(format!("{bmp_base}{sfx}"), idx as i32);
        }
        for (idx, sfx) in suffixes.iter().enumerate() {
            m.insert(format!("{supp_base}{sfx}"), 100 + idx as i32);
        }
        assert!(m.treeified());
        assert_eq!(m.capacity(), 64);
        for (idx, sfx) in suffixes.iter().enumerate() {
            assert_eq!(m.get(&format!("{bmp_base}{sfx}")), Some(&(idx as i32)));
            assert_eq!(
                m.get(&format!("{supp_base}{sfx}")),
                Some(&(100 + idx as i32))
            );
        }

        let actual_tags: Vec<String> = m
            .keys()
            .map(|k| {
                if let Some(rest) = k.strip_prefix(bmp_base) {
                    format!("B:{rest}")
                } else if let Some(rest) = k.strip_prefix(supp_base) {
                    format!("S:{rest}")
                } else {
                    k.clone()
                }
            })
            .collect();
        assert_eq!(
            actual_tags,
            vec![
                "B:AaAaBB", "B:AaAaAa", "B:AaBBAa", "B:AaBBBB", "B:BBAaAa", "B:BBAaBB", "S:AaAaAa",
                "S:AaAaBB", "S:AaBBAa", "S:AaBBBB", "S:BBAaAa", "S:BBAaBB",
            ]
        );
    }

    fn max_tree_depth<K, V>(map: &JavaHashMap<K, V>) -> usize
    where
        K: Clone + JavaHash + Eq + PartialEq + Hash + Debug,
        V: Clone + Debug + Default,
    {
        fn subtree_depth(links: &[TreeLink], node: u32) -> usize {
            if node == NONE {
                return 0;
            }
            let lk = links[node as usize];
            1 + subtree_depth(links, lk.left).max(subtree_depth(links, lk.right))
        }
        let Some(links) = map.tree_links.as_deref() else {
            return 0;
        };
        map.buckets
            .iter()
            .copied()
            .filter(|&h| h != NONE && map.is_tree_node(h))
            .map(|h| subtree_depth(links, h))
            .max()
            .unwrap_or(0)
    }

    fn count_lookup_probe_steps<K, V>(map: &JavaHashMap<K, V>, key: &K) -> (bool, usize)
    where
        K: Clone + JavaHash + Eq + PartialEq + Hash + Debug,
        V: Clone + Debug + Default,
    {
        if map.buckets.is_empty() {
            return (false, 0);
        }
        let hash = JavaHashMap::<K, V>::spread(key);
        let idx = (hash as usize) & (map.buckets.len() - 1);
        let first = map.buckets[idx];
        if first == NONE {
            return (false, 0);
        }
        let mut steps = 1usize;
        if &map.nodes[first as usize].key == key {
            return (true, steps);
        }
        let mut curr = map.nodes[first as usize].next;
        if curr == NONE {
            return (false, steps);
        }
        if map.is_tree_node(first) {
            let links = map.tree_links.as_deref().expect("tree_links present");
            let mut p = if links[first as usize].parent != NONE {
                map.tree_root(first)
            } else {
                first
            };
            let h = hash as i32;
            while p != NONE {
                steps += 1;
                let node = &map.nodes[p as usize];
                let link = links[p as usize];
                let ph = JavaHashMap::<K, V>::spread(&node.key) as i32;
                if ph > h {
                    p = link.left;
                } else if ph < h {
                    p = link.right;
                } else if &node.key == key {
                    return (true, steps);
                } else if link.left == NONE {
                    p = link.right;
                } else if link.right == NONE {
                    p = link.left;
                } else {
                    p = match key.java_compare(&node.key) {
                        std::cmp::Ordering::Less => link.left,
                        std::cmp::Ordering::Equal => return (true, steps),
                        std::cmp::Ordering::Greater => link.right,
                    };
                }
            }
            return (false, steps);
        }
        while curr != NONE {
            steps += 1;
            if &map.nodes[curr as usize].key == key {
                return (true, steps);
            }
            curr = map.nodes[curr as usize].next;
        }
        (false, steps)
    }

    #[test]
    fn ten_thousand_stride_collisions_run_in_n_log_n_with_openjdk_order() {
        // Deterministic O(log N) / O(N log N) complexity verification across N in {2_500, 5_000, 10_000}
        // (no host-dependent wall-clock gate per AGENTS.md §14).
        let mut avg_steps_by_scale = Vec::new();
        let mut final_map: Option<JavaHashMap<i32, i32>> = None;

        for &n in &[2_500usize, 5_000, 10_000] {
            let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
            for i in 0..(n as i32) {
                // spread(k) = (i << 16) | 0x1234 -> all keys collide into bucket 0x1234 at cap <= 65,536!
                let k = (i << 16) | ((i ^ 0x1234) & 0xFFFF);
                m.insert(k, i);
            }
            assert!(m.treeified());
            assert_eq!(m.len(), n);
            assert_openjdk_invariants(&m);

            let ceil_log2 = (usize::BITS - n.leading_zeros()) as usize;
            let rb_depth_bound = 2 * ceil_log2;
            let depth = max_tree_depth(&m);
            assert!(
                depth <= rb_depth_bound,
                "Red-Black depth {depth} exceeded 2*ceil(log2({n}+1))={rb_depth_bound}"
            );

            let mut total_steps = 0usize;
            for i in 0..(n as i32) {
                let k = (i << 16) | ((i ^ 0x1234) & 0xFFFF);
                assert_eq!(m.get(&k), Some(&i));
                let (found, steps) = count_lookup_probe_steps(&m, &k);
                assert!(found);
                assert!(
                    steps <= rb_depth_bound + 1,
                    "single lookup took {steps} steps at N={n} (bound {})",
                    rb_depth_bound + 1
                );
                total_steps += steps;
            }
            let total_bound = n * rb_depth_bound;
            assert!(
                total_steps <= total_bound,
                "total lookup probes {total_steps} exceeded O(N log N) bound {total_bound} at N={n}"
            );
            avg_steps_by_scale.push((total_steps as f64) / (n as f64));
            if n == 10_000 {
                final_map = Some(m);
            }
        }

        // Verify average probe steps grow logarithmically (~+1 step per 2x doubling, strictly <= +3.0):
        for pair in avg_steps_by_scale.windows(2) {
            let delta = pair[1] - pair[0];
            assert!(
                delta > 0.0 && delta <= 3.0,
                "expected logarithmic step growth per 2x N doubling, got {avg_steps_by_scale:?}"
            );
        }

        let mut m = final_map.expect("10,000-key map constructed");
        let first_five_before: Vec<i32> = m.keys().take(5).copied().collect();
        assert_eq!(
            first_five_before,
            vec![134_157_771, 67_047_883, 33_493_963, 16_716_491, 8_327_755]
        );
        for i in (0..10_000i32).step_by(2) {
            let k = (i << 16) | ((i ^ 0x1234) & 0xFFFF);
            assert_eq!(m.remove(&k), Some(i));
        }
        assert_eq!(m.len(), 5_000);
        assert!(m.treeified());
        assert_openjdk_invariants(&m);
        assert!(max_tree_depth(&m) <= 26);
        let first_five_after: Vec<i32> = m.keys().take(5).copied().collect();
        assert_eq!(
            first_five_after,
            vec![268_377_547, 134_157_771, 67_047_883, 33_493_963, 16_716_491]
        );
    }

    #[test]
    fn iter_clone_without_clone_value_and_duplicate_key_order() {
        #[derive(Debug, PartialEq, Eq)]
        struct NonCloneVal(i32);

        let nodes = [
            Node {
                key: 10,
                val: NonCloneVal(100),
                next: 2,
            },
            Node {
                key: 20,
                val: NonCloneVal(200),
                next: NONE,
            },
            Node {
                key: 26,
                val: NonCloneVal(260),
                next: NONE,
            },
        ];
        let buckets = [NONE, 0, 1, NONE];
        let mut iter: Iter<'_, i32, NonCloneVal> = Iter {
            buckets: &buckets,
            nodes: &nodes,
            bucket_idx: 0,
            curr: NONE,
        };
        assert_eq!(iter.next(), Some((&10, &NonCloneVal(100))));
        let cloned = iter.clone();
        let rest_orig: Vec<(i32, i32)> = iter.map(|(&k, v)| (k, v.0)).collect();
        let rest_cloned: Vec<(i32, i32)> = cloned.map(|(&k, v)| (k, v.0)).collect();
        assert_eq!(rest_orig, vec![(26, 260), (20, 200)]);
        assert_eq!(rest_orig, rest_cloned);

        // Verify duplicate keys retain first-seen reference identity and unsized `str` works:
        let s1 = String::from("chr1");
        let s2 = String::from("chrX");
        let s1_dup = String::from("chr1");
        let input = [&s1, &s2, &s1_dup];
        let ordered = java_hashmap_key_order(input.iter().copied());
        assert_eq!(ordered.len(), 2);
        assert!(std::ptr::eq(ordered[0], &s2));
        assert!(std::ptr::eq(ordered[1], &s1));

        let str_refs: [&str; 4] = ["chr1", "chrX", "chr1", "chrX"];
        assert_eq!(java_hashmap_key_order(str_refs), vec!["chrX", "chr1"]);
    }

    #[test]
    fn test_javamap_size_into_entries_and_ascii_hash_order() {
        assert!(
            std::mem::size_of::<JavaHashMap<i32, crate::variations::Vars>>() <= 96,
            "expected JavaHashMap<i32, Vars> <= 96 bytes, got {}",
            std::mem::size_of::<JavaHashMap<i32, crate::variations::Vars>>()
        );

        let samples = [
            "",
            "A",
            "+ACGT",
            "-12",
            "chr1",
            "chr10",
            "chrX",
            "\u{00E9}",
            "A\u{FFFD}",
            "\u{E000}",
            "\u{FFFD}",
            "\u{10000}",
            "\u{1F600}",
        ];
        for s in samples {
            let expected_hash = s
                .encode_utf16()
                .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(i32::from(c)));
            assert_eq!(
                java_string_hash(s),
                expected_hash,
                "hash mismatch for {s:?}"
            );
            for t in samples {
                let expected_cmp = s.encode_utf16().cmp(t.encode_utf16());
                assert_eq!(
                    s.java_compare(t),
                    expected_cmp,
                    "java_compare mismatch for {s:?} vs {t:?}"
                );
            }
        }

        // Empty map
        let empty: JavaHashMap<i32, i32> = JavaHashMap::new();
        assert_eq!(
            empty.clone().into_entries().collect::<Vec<_>>(),
            empty.into_iter().collect::<Vec<_>>()
        );

        // Linked map (with removals leaving free slots and trailing empty buckets)
        let mut linked: JavaHashMap<i32, i32> = JavaHashMap::new();
        for k in [3, 19, 35, 7, 23, 100, 200, 300] {
            linked.insert(k, k * 10);
        }
        assert_eq!(linked.remove(&19), Some(190));
        assert!(!linked.treeified());
        assert_eq!(
            linked.clone().into_entries().collect::<Vec<_>>(),
            linked.into_iter().collect::<Vec<_>>()
        );

        // Treeified (2^16 + 1)-stride map (cancels `h ^ (h >>> 16)` in low 16 bits -> bucket 0)
        let mut tree: JavaHashMap<i32, i32> = JavaHashMap::new();
        for i in 0..20i32 {
            tree.insert((i << 16) | i, i * 100);
        }
        assert!(tree.treeified());
        assert_eq!(tree.remove(&((3 << 16) | 3)), Some(300));
        assert_eq!(
            tree.clone().into_entries().collect::<Vec<_>>(),
            tree.into_iter().collect::<Vec<_>>()
        );
    }
}
