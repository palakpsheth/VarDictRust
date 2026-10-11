use std::fmt::Debug;
use std::hash::Hash;

const DEFAULT_CAPACITY: usize = 16;
const TREEIFY_THRESHOLD: usize = 8;
const MIN_TREEIFY_CAPACITY: usize = 64;
const DEFAULT_LOAD_FACTOR: f32 = 0.75;

#[derive(Clone, Debug)]
pub struct JavaHashMap<K, V> {
    table: Vec<Vec<(K, V)>>,
    size: usize,
    threshold: usize,
    treeified: bool,
}

impl<K, V> JavaHashMap<K, V>
where
    K: Copy + Into<i32> + Eq + PartialEq + Hash + Debug,
    V: Clone + Debug + Default,
{
    pub fn new() -> Self {
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

    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
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

    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
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

    pub fn compute_if_absent<F>(&mut self, key: K, f: F) -> &mut V
    where
        F: FnOnce(&K) -> V,
    {
        if self.table.is_empty() || self.size > self.threshold {
            // OpenJDK HashMap.computeIfAbsent pre-resize check (size > threshold)
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
                // OpenJDK HashMap.computeIfAbsent inserts at head of bin
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

        // Determine target bin and position via immutable borrow
        let pos = self
            .table
            .get(idx)
            .and_then(|bin| bin.iter().position(|(k, _)| *k == key));

        let (target_bin, target_pos) = match pos {
            Some(p) => (idx, p),
            None => {
                // If not found at idx (e.g. treeify_bin resized), scan table
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

        // Fallback for compiler completeness
        if self.table.is_empty() {
            self.table.push(Vec::new());
        }
        let bin = &mut self.table[0];
        if bin.is_empty() {
            bin.push((key, V::default()));
        }
        &mut bin[0].1
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
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

    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.table.iter().flat_map(|bin| bin.iter().map(|(k, _)| k))
    }

    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.table
            .iter()
            .flat_map(|bin| bin.iter().map(|(k, v)| (k, v)))
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

    fn into_iter(self) -> Self::IntoIter {
        let mut items = Vec::new();
        for mut bin in self.table {
            items.append(&mut bin);
        }
        items.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
