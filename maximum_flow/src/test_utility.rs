use crate::graph::Flow;
use std::fs::read_to_string;
use std::path::PathBuf;

pub struct GraphInstance {
    pub num_nodes: usize,
    pub edges: Vec<(usize, usize, Flow)>,
    pub source: usize,
    pub sink: usize,
}

pub fn read_expected(file_path: &PathBuf) -> Flow {
    let data = read_to_string(file_path).unwrap();
    data.trim().parse().unwrap()
}

pub fn read_instance_aoj(file_path: &PathBuf) -> GraphInstance {
    let data = read_to_string(file_path).unwrap();
    let data: Vec<&str> = data.trim().split('\n').collect();

    let (num_nodes, source, sink) = {
        let header: Vec<&str> = data[0].trim().split(' ').collect();
        let num_nodes = header[0].parse().unwrap();
        (num_nodes, 0, num_nodes - 1)
    };

    let mut edges = Vec::new();
    for line in data[1..].iter() {
        let line: Vec<&str> = line.trim().split(' ').collect();
        let (from, to, capacity) = (
            line[0].parse().unwrap(),
            line[1].parse().unwrap(),
            line[2].parse().unwrap(),
        );
        edges.push((from, to, capacity));
    }

    GraphInstance {
        num_nodes,
        edges,
        source,
        sink,
    }
}

pub fn read_instance_libreoj(file_path: &PathBuf) -> GraphInstance {
    let data = read_to_string(file_path).unwrap();
    let data: Vec<&str> = data.trim().split('\n').collect();

    let (num_nodes, source, sink) = {
        let header: Vec<&str> = data[0].trim().split(' ').collect();
        let (num_nodes, source, sink): (usize, usize, usize) = (
            header[0].parse().unwrap(),
            header[2].parse().unwrap(),
            header[3].parse().unwrap(),
        );
        (num_nodes, source - 1, sink - 1)
    };

    let mut edges = Vec::new();
    for line in data[1..].iter() {
        let line: Vec<&str> = line.trim().split(' ').collect();
        let (from, to, capacity): (usize, usize, Flow) = (
            line[0].parse().unwrap(),
            line[1].parse().unwrap(),
            line[2].parse().unwrap(),
        );
        edges.push((from - 1, to - 1, capacity));
    }

    GraphInstance {
        num_nodes,
        edges,
        source,
        sink,
    }
}
