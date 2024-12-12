use std::collections::HashMap;

/// Implementation of a Trie, where all items in the Trie *must* be of the same length
struct TrieNode {
    children: HashMap<u8, Box<TrieNode>> // need to use boxes since recursive so don't know size. HashMap has known size (it's a smart pointer)
}

impl TrieNode {

    pub fn new() -> TrieNode {
        TrieNode {children: HashMap::new()}
    }

    pub fn is_end_of_word(&self) -> bool {
        // since all items in the Trie must be of the same length, this can be inferred based on the presence/absence of children
        self.children.is_empty()
    }

    pub fn has_child(&self, byte: u8) -> bool {
        self.children.contains_key(&byte)
    }
    
    pub fn add_child(&mut self, byte: u8) {
        let n = TrieNode::new();
        self.children.insert(byte, Box::new(n));
    }

    pub fn get_child(&self, byte: u8) -> Option<&Box<TrieNode>> {
        self.children.get(&byte)
    }

    pub fn get_child_mut(&mut self, byte: u8) -> Option<&mut Box<TrieNode>> {
        self.children.get_mut(&byte)
    }

    pub fn get_children_ids(&self) -> Vec<u8> {
        let v: Vec<u8> = self.children.keys().cloned().collect();
        v
    }

}

pub struct Trie {
    word_length: usize,
    word_count: usize,
    root: TrieNode
}

impl Trie {
    
    pub fn new () -> Trie {
        Trie {root: TrieNode::new(), word_count: 0, word_length: 0}
    }
    
    pub fn contains_word(&self, word: &[u8]) -> bool {
        //! Check if the Trie contains a given word.
        //! # Examples
        //! ```
        //! use barcodes::trie::Trie;
        //! let mut t = Trie::new();
        //! t.add_word(b"hello");
        //! assert_eq!(t.contains_word(b"hello"), true);
        //! assert_eq!(t.contains_word(b"goodbye"), false);
        //! ```
        let mut node = &self.root;

        for &byte in word {
            if !node.has_child(byte) {
                return false;
            } else {
                node = node.get_child(byte).unwrap();
            }
        }

        true
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let mut node = &mut (self.root);
        let mut word_is_new = false; // flipped to true once a new child is added. Cannot just use self.contains_word(word) because that requires an immutable ref, and we already have a mutable one

        if self.word_count != 0 {
            assert_eq!(self.word_length, word.len());
        } else {
            self.word_length = word.len();
        }

        for &byte in word {
            if !node.has_child(byte) {
                node.add_child(byte);
                word_is_new = true;
            }
            node = node.get_child_mut(byte).unwrap();
        }

        self.word_count += if word_is_new {1} else {0};
    }

    /// Get the number of items in the Trie
    pub fn len(&self) -> usize {
        self.word_count
    }

    /// Get the length of an item in the Trie (all items are of the same length)
    pub fn word_length(&self) -> usize {
        if self.word_count == 0 {
            panic!("word_length is undefined on an empty Trie");
        }

        self.word_length
    }

    pub fn get_words_within_hamming_distance(&self, word: &[u8], max_distance: usize) -> Vec<(String, usize)> {
        assert_eq!(word.len(), self.word_length);
        self._get_within_hamming_distance(&self.root, word, &String::from(""), 0, max_distance)
    }

    fn _get_within_hamming_distance(&self, node: &TrieNode, word: &[u8], prefix: &String, current_distance: usize, max_distance: usize) -> Vec<(String, usize)> {
        let mut matches: Vec<(String, usize)> = Vec::new();
        
        if node.is_end_of_word() {
            assert!(current_distance <= max_distance);
            matches.push((prefix.clone(), current_distance));
        } else {
            for child_id in node.get_children_ids() {
                let cost = if child_id == word[0] {0} else {1};
                if current_distance + cost > max_distance {
                    continue
                }
                let child_node = node.get_child(child_id).unwrap();
                let mut new_prefix = prefix.clone();
                new_prefix.push(child_id as char);
                let mut m = self._get_within_hamming_distance(child_node, &word[1..word.len()], &new_prefix, current_distance + cost, max_distance);
                matches.append(&mut m);
            }
        }

        matches
    }


}