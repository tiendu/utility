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
}

const M: isize = 2;
const MM: isize = -1;
const G: isize = -2;
const E: isize = -1;

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let input = args[1].to_string();
    let file = File::open(&input.to_string()).unwrap();
    let reader = BufReader::new(file);
    let mut id = String::new();
    let mut storage: Vec<IdSequence> = Vec::new();
    for line in reader.lines() {
        let l = &mut line.as_ref().unwrap().to_string();
        let s = l.split("");
        let vec: Vec<&str> = s.collect();
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
    for i in storage {
        sequences.push(i.sequence);
    }
    let mut sequences: Vec<&str> = sequences
        .iter()
        .map(|s| &**s)
        .collect();
    sequences.sort_by(|a, b| a.len().cmp(&b.len()));
    sequences.reverse();
    let result = allpairs_cmp(sequences);
    println!("{:?}", result);
    Ok(())
}

fn create_reverse_complement(seq: &str) -> String {
    let mut temp = String::from("");
    let lst: Vec<&str> = seq
        .split("")
        .collect();
    let mut count = lst.len() - 2;
    while count > 0 {
        if lst[count] == "A" {
            temp += &"T".to_owned();
        } else if lst[count] == "T" {
            temp += &"A".to_owned();
        } else if lst[count] == "G" {
            temp += &"C".to_owned();
        } else if lst[count] == "C" {
            temp += &"G".to_owned();
        }
        count -= 1;
    }
    return temp
}

fn align(seq1: &str, seq2: &str, bound: isize) -> ([String; 2], [usize; 2], [usize; 2], [isize; 3]) {
    let len1 = seq1.len();
    let len2 = seq2.len();
    let mut mat: Vec<Vec<(isize, isize)>> =  vec![vec![(bound, 0); len2 + 1]; len1 + 1];
    let mut lst1: Vec<&str> = seq1
        .split("")
        .collect();
    let mut lst2: Vec<&str> = seq2
        .split("")
        .collect();
    lst1 = lst1[1 .. ].to_vec();
    lst2 = lst2[1 .. ].to_vec();
    let mut max_score = (0, 0);
    let mut end_pos = [0, 0];
    let left: isize = 1;
    let diag: isize = 2;
    let up: isize = 3;
    for i in 1 .. len1 + 1 {
        for j in 1 .. len2 + 1 {
            let p = if lst1[i - 1] == lst2[j - 1] { M } else { MM };
            let score = *[
                (0, 0),
                (mat[i - 1][j].0 + G, up),
                (mat[i][j - 1].0 + G, left),
                (mat[i - 1][j - 1].0 + p, diag),
            ]
                .iter()
                .max_by(|x, y| x.0.cmp(&y.0))
                .unwrap();
            if score.0 >= max_score.0 {
                max_score = score;
                end_pos = [i, j];
            }
            if score.0 + M * (*[len1 - i, len2 - j].iter().max().unwrap() as isize) >= bound {
                mat[i][j] = score;   
            }
        }
    }
    let mut x = end_pos[0];
    let mut y = end_pos[1];
    let mut aln1: Vec<&str> = Vec::new();
    let mut aln2: Vec<&str> = Vec::new();
    let mut start_pos = [0, 0];
    let mut _match: isize = 0;
    let mut _mismatch: isize = 0;
    let mut _gap: isize = 0;
    loop {
        if mat[x][y].1 == 1 && mat[x][y].0 != bound { 
            y -= 1;
            aln1.push("-");
            aln2.push(lst2[y]);
            _gap += 1;
        } else if mat[x][y].1 == 2 && mat[x][y].0 != bound {
            x -= 1;
            y -= 1;
            aln1.push(lst1[x]);
            aln2.push(lst2[y]);
            if lst1[x] == lst2[y] {
                _match += 1;
            } else {
                _mismatch += 1;
            }
        } else if mat[x][y].1 == 3 && mat[x][y].0 != bound {
            x -= 1;
            aln1.push(lst1[x]);
            aln2.push("-");
            _gap += 1;
        } else if mat[x][y].0 == bound {
            start_pos = [x, y];
            break;
        }
    }
    aln1.reverse();
    aln2.reverse();
    let score = [_match, _mismatch, _gap];
    return ([aln1.join("").to_string(), aln2.join("").to_string()], [start_pos[0] + 1, end_pos[0]], [start_pos[1] + 1, end_pos[1]], score);
}

fn compute_bound(seq1: &str, seq2: &str, seq3: &str) -> isize {
    let mut bound: isize = 0;
    let res1 = align(seq1, seq2, bound);
    let res2 = align(seq1, seq3, bound);
    if res1.1[1] >= res2.1[0] && res2.1[1] >= res1.1[0] {
        let _mismatch: isize = MM * ((res1.3[1] + res2.3[1]) as isize);
        let _gap: isize = *[G +  E * ((res1.3[2] + res2.3[2] - 1) as isize), G * ((res1.3[2] + res2.3[2]) as isize)]
            .iter()
            .max()
            .unwrap();
        let _match: isize = (*[res1.1[1], res2.1[1]]
            .iter()
            .min()
            .unwrap() as isize) - (*[res1.1[0], res2.1[0]]
                .iter()
                .max()
                .unwrap() as isize) + 1 - ((res1.3[1] + res2.3[1]) as isize);
        bound = *[0, M * _match + _mismatch + _gap]
            .iter()
            .max()
            .unwrap();
    }
    return bound;
}

fn allpairs_cmp(seqs: Vec<&str>) -> Vec<([String; 2], [usize; 2], [usize; 2], [isize; 3])> {
    let n = seqs.len();
    let mut res: Vec<([String; 2], [usize; 2], [usize; 2], [isize; 3])> = Vec::new();
    for a in 1 .. n - 1 {
        for b in a + 1 .. n {
            let mut bound: isize = 0;
            for c in 1 .. a - 1 {
                bound = *[bound, compute_bound(seqs[c], seqs[a], seqs[b])]
                    .iter()
                    .max()
                    .unwrap();
            }
            res.push(align(seqs[a], seqs[b], bound));
        }
    }
    return res;
}
