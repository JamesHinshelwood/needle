use std::cmp::max;
use std::fmt;

use nalgebra::DMatrix;

fn needleman_wunsch<T, F>(first: Vec<T>, second: Vec<T>, similarity_fn: F, indel_score: i32) -> (Vec<Option<T>>, Vec<Option<T>>)
where F: Fn(T, T) -> i32,
      T: Copy {
    let mut matrix: DMatrix<i32> = DMatrix::zeros(first.len() + 1, second.len() + 1);

    for (i, element) in matrix.row_mut(0).iter_mut().enumerate() {
        *element = indel_score * (i as i32);
    }

    for (i, element) in matrix.column_mut(0).iter_mut().enumerate() {
        *element = indel_score * (i as i32);
    }

    for i in 1..matrix.nrows() {
        for j in 1..matrix.ncols() {
            let mismatch = matrix[(i - 1, j - 1)] + similarity_fn(first[i - 1], second[j - 1]);
            let delete = matrix[(i - 1, j)] + indel_score;
            let insert = matrix[(i, j - 1)] + indel_score;

            matrix[(i, j)] = max(max(mismatch, delete), insert);
        }
    }

    let mut alignment_a = Vec::with_capacity(max(first.len(), second.len()));
    let mut alignment_b = Vec::with_capacity(max(first.len(), second.len()));

    let mut i = matrix.nrows() - 1;
    let mut j = matrix.ncols() - 1;
    while i > 0 || j > 0 {
        if i > 0 && j > 0 && matrix[(i, j)] == matrix[(i - 1, j - 1)] + similarity_fn(first[i - 1], second[j - 1]) {
            alignment_a.push(Some(first[i - 1]));
            alignment_b.push(Some(second[j - 1]));
            i = i - 1;
            j = j - 1;
        }
        else if i > 0 && matrix[(i, j)] == matrix[(i - 1, j)] + indel_score {
            alignment_a.push(Some(first[i - 1]));
            alignment_b.push(None);
            i = i - 1;
        }
        else {
            alignment_a.push(None);
            alignment_b.push(Some(second[j - 1]));
            j = j - 1;
        }
    }

    alignment_a.reverse();
    alignment_b.reverse();

    (alignment_a, alignment_b)
}

#[allow(dead_code)]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum AminoAcid {
        A,
        R,
        N,
        D,
        C,
        Q,
        E,
        G,
        H,
        I,
        L,
        K,
        M,
        F,
        P,
        S,
        T,
        W,
        Y,
        V,
        B,
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

fn blosum_50(a: AminoAcid, b: AminoAcid) -> i32 {
    use self::AminoAcid::*;
    match a {
        A => { match b { A => 5, R => -2, N => -1, D => -2, C => -1, Q => -1, E => -1, G => 0, H => -2, I => -1, L => -2, K => -1, M => -1, F => -3, P => -1, S => 1, T => 0, W => -3, Y => -2, V => 0, B => -2, }},
        R => { match b { A => -2, R => 7, N => -1, D => -2, C => -4, Q => 1, E => 0, G => -3, H => 0, I => -4, L => -3, K => 3, M => -2, F => -3, P => -3, S => -1, T => -1, W => -3, Y => -1, V => -3, B => -1, }},
        N => { match b { A => -1, R => -1, N => 7, D => 2, C => -2, Q => 0, E => 0, G => 0, H => 1, I => -3, L => -4, K => 0, M => -2, F => -4, P => -2, S => 1, T => 0, W => -4, Y => -2, V => -3, B => 4, }},
        D => { match b { A => -2, R => -2, N => 2, D => 8, C => -4, Q => 0, E => 2, G => -1, H => -1, I => -4, L => -4, K => -1, M => -4, F => -5, P => -1, S => 0, T => -1, W => -5, Y => -3, V => -4, B => 5, }},
        C => { match b { A => -1, R => -4, N => -2, D => -4, C => 13, Q => -3, E => -3, G => -3, H => -3, I => -2, L => -2, K => -3, M => -2, F => -2, P => -4, S => -1, T => -1, W => -5, Y => -3, V => -1, B => -3, }},
        Q => { match b { A => -1, R => 1, N => 0, D => 0, C => -3, Q => 7, E => 2, G => -2, H => 1, I => -3, L => -2, K => 2, M => 0, F => -4, P => -1, S => 0, T => -1, W => -1, Y => -1, V => -3, B => 0, }},
        E => { match b { A => -1, R => 0, N => 0, D => 2, C => -3, Q => 2, E => 6, G => -3, H => 0, I => -4, L => -3, K => 1, M => -2, F => -3, P => -1, S => -1, T => -1, W => -3, Y => -2, V => -3, B => 1, }},
        G => { match b { A => 0, R => -3, N => 0, D => -1, C => -3, Q => -2, E => -3, G => 8, H => -2, I => -4, L => -4, K => -2, M => -3, F => -4, P => -2, S => 0, T => -2, W => -3, Y => -3, V => -4, B => -1, }},
        H => { match b { A => -2, R => 0, N => 1, D => -1, C => -3, Q => 1, E => 0, G => -2, H => 10, I => -4, L => -3, K => 0, M => -1, F => -1, P => -2, S => -1, T => -2, W => -3, Y => 2, V => -4, B => 0, }},
        I => { match b { A => -1, R => -4, N => -3, D => -4, C => -2, Q => -3, E => -4, G => -4, H => -4, I => 5, L => 2, K => -3, M => 2, F => 0, P => -3, S => -3, T => -1, W => -3, Y => -1, V => 4, B => -4, }},
        L => { match b { A => -2, R => -3, N => -4, D => -4, C => -2, Q => -2, E => -3, G => -4, H => -3, I => 2, L => 5, K => -3, M => 3, F => 1, P => -4, S => -3, T => -1, W => -2, Y => -1, V => 1, B => -4, }},
        K => { match b { A => -1, R => 3, N => 0, D => -1, C => -3, Q => 2, E => 1, G => -2, H => 0, I => -3, L => -3, K => 6, M => -2, F => -4, P => -1, S => 0, T => -1, W => -3, Y => -2, V => -3, B => 0, }},
        M => { match b { A => -1, R => -2, N => -2, D => -4, C => -2, Q => 0, E => -2, G => -3, H => -1, I => 2, L => 3, K => -2, M => 7, F => 0, P => -3, S => -2, T => -1, W => -1, Y => 0, V => 1, B => -3, }},
        F => { match b { A => -3, R => -3, N => -4, D => -5, C => -2, Q => -4, E => -3, G => -4, H => -1, I => 0, L => 1, K => -4, M => 0, F => 8, P => -4, S => -3, T => -2, W => 1, Y => 4, V => -1, B => -4, }},
        P => { match b { A => -1, R => -3, N => -2, D => -1, C => -4, Q => -1, E => -1, G => -2, H => -2, I => -3, L => -4, K => -1, M => -3, F => -4, P => 10, S => -1, T => -1, W => -4, Y => -3, V => -3, B => -2, }},
        S => { match b { A => 1, R => -1, N => 1, D => 0, C => -1, Q => 0, E => -1, G => 0, H => -1, I => -3, L => -3, K => 0, M => -2, F => -3, P => -1, S => 5, T => 2, W => -4, Y => -2, V => -2, B => 0, }},
        T => { match b { A => 0, R => -1, N => 0, D => -1, C => -1, Q => -1, E => -1, G => -2, H => -2, I => -1, L => -1, K => -1, M => -1, F => -2, P => -1, S => 2, T => 5, W => -3, Y => -2, V => 0, B => 0, }},
        W => { match b { A => -3, R => -3, N => -4, D => -5, C => -5, Q => -1, E => -3, G => -3, H => -3, I => -3, L => -2, K => -3, M => -1, F => 1, P => -4, S => -4, T => -3, W => 15, Y => 2, V => -3, B => -5, }},
        Y => { match b { A => -2, R => -1, N => -2, D => -3, C => -3, Q => -1, E => -2, G => -3, H => 2, I => -1, L => -1, K => -2, M => 0, F => 4, P => -3, S => -2, T => -2, W => 2, Y => 8, V => -1, B => -3, }},
        V => { match b { A => 0, R => -3, N => -3, D => -4, C => -1, Q => -3, E => -3, G => -4, H => -4, I => 4, L => 1, K => -3, M => 1, F => -1, P => -3, S => -2, T => 0, W => -3, Y => -1, V => 5, B => -4, }},
        B => { match b { A => -2, R => -1, N => 4, D => 5, C => -3, Q => 0, E => 1, G => -1, H => 0, I => -4, L => -4, K => 0, M => -3, F => -4, P => -2, S => 0, T => 0, W => -5, Y => -3, V => -4, B => 5, }},
    }
}

fn main() {
    use self::AminoAcid::*;
    let first  = vec![P, A, W, H, E, A, E];
    let second = vec![H, E, A, G, A, W, G, H, E, E];

    let (a, b) = needleman_wunsch(first, second, blosum_50, -8);

    println!("{:?}", a.iter().map(|c| c.map_or('-', |x| x.to_string().chars().next().unwrap())).collect::<String>());
    println!("{:?}", b.iter().map(|c| c.map_or('-', |x| x.to_string().chars().next().unwrap())).collect::<String>());
}