use std::env;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

struct IdSequence {
    id: String,
    sequence: String,
}

impl IdSequence {
    fn new(id: String, sequence: String) -> IdSequence {
        IdSequence { id: id, sequence: sequence }
    }
    fn base_frequency(&self) -> Vec<f64> {
        let bases: Vec<String> = vec!["A", "T", "G", "C"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut base_freq: Vec<f64> = Vec::new();
        for i in bases {
            base_freq.push((self.sequence
                .matches(&i)
                .count() as f64) / (self.sequence
                    .chars()
                    .count()) as f64);
        }
        return base_freq;
    }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let input = args[1].to_string();
    let file = File::open(&input.to_string()).unwrap();
    let reader = BufReader::new(file);
    let mut id = String::new();
    let mut storage: Vec<IdSequence> = Vec::new();
    for line in reader.lines() {
        let l = &mut line
            .as_ref()
            .unwrap()
            .to_string();
        let vec: Vec<&str> = l
            .split("")
            .collect();
        let mut name = String::new();
        if vec[1] == ">" {
            for i in 1 .. vec.len() {
                if vec[i] == " " {
                    name = vec[2 .. i].join("");
                    break;
                } else if vec[i] != "" { 
                    name = vec[2 .. vec.len()].join("");
                }
            }
            id = name;
        } else {
            l.retain(|c| !c.is_whitespace());
            storage.push(IdSequence::new(id.to_string(), l.to_string()));
        }
    }
    let mut sequences: Vec<String> = Vec::new();
    let mut gc_content: f64 = 0.0;
    for i in &storage {
        sequences.push(i.sequence.clone());
        gc_content += i.base_frequency()[2] + i.base_frequency()[3];
    }
    gc_content = gc_content * 100.0 / (sequences.len() as f64);
    let n50 = Nx(sequences.clone(), 50);
    let n90 = Nx(sequences.clone(), 90);
    let totsum: usize = sequences
        .clone()
        .iter()
        .map(|x| x.len())
        .sum();
    let sqsum: usize = sequences
        .clone()
        .iter()
        .map(|x| x.len().pow(2))
        .sum();
    let auN: f64 = (sqsum as f64) / (totsum as f64);
    println!("{}\tN50: {}\tN90: {}\tauN: {:.3}\tGC: {:.3}", input, n50, n90, auN, gc_content);
    Ok(())
}

fn Nx(mut seqs: Vec<String>, x: usize) -> usize {
    seqs.sort_by(|a, b| a.len().cmp(&b.len()));
    let mut cumsum: usize = 0;
    let totsum: usize = seqs
        .iter()
        .map(|x| x.len())
        .sum();
    for i in 0 .. seqs.len() {
        cumsum += seqs[i].len();
        if 100 * cumsum >= totsum * (100 - x) {
            return seqs[i].len()
        }
    }
    return 0
} 
