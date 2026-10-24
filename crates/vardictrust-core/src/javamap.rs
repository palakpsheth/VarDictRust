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

#[derive(Clone, Debug)]
pub struct JavaHashMap<K, V> {
    buckets: Vec<u32>,
    nodes: Vec<Node<K, V>>,
    free_head: u32,
    size: usize,
    threshold: usize,
    treeified: bool,
    fallback: V,
}

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

impl<K, V> JavaHashMap<K, V>
where
    K: Copy + Into<i32> + Eq + PartialEq + Hash + Debug,
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
        }
    }

    fn spread(key: K) -> u32 {
        let h = key.into() as u32;
        h ^ (h >> 16)
    }

    fn alloc_node(&mut self, key: K, val: V, next: u32) -> u32 {
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

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    pub fn contains_key(&self, key: &K) -> bool {
        self.get(key).is_some()
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

        if bin_len + 1 >= TREEIFY_THRESHOLD {
            self.treeify_bin();
        }
        self.size += 1;
        if self.size > self.threshold {
            self.resize();
        }
        None
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

    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
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
    K: Copy + Into<i32> + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<K, V> IntoIterator for JavaHashMap<K, V>
where
    K: Copy + Into<i32> + Eq + PartialEq + Hash + Debug,
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
                if bin_len >= TREEIFY_THRESHOLD {
                    self.treeify_bin();
                }
                self.size += 1;
                if self.size > self.threshold {
                    self.resize();
                }
            }
            None
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
}
