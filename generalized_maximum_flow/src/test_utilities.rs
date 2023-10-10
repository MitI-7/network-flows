use crate::graph::Flow;
use std::fs::read_to_string;
use std::path::{Path, PathBuf};

pub struct GraphInstancve {
    pub num_nodes: usize,
    pub edges: Vec<(usize, usize, Flow, Flow)>,
    pub source: usize,
    pub sink: usize,
}

pub fn read_expected(file_path: &Path) -> Flow {
    let data = read_to_string(file_path).unwrap();
    data.trim().parse().unwrap()
}

pub fn read_graph_instance(file_path: &PathBuf) -> GraphInstancve {
    let data = read_to_string(file_path).unwrap();
    let data: Vec<&str> = data.trim().split('\n').collect();

    let (num_nodes, source, sink) = {
        let header: Vec<&str> = data[0].trim().split(' ').collect();
        let num_nodes = header[0].parse().unwrap();
        let source = header[2].parse().unwrap();
        let sink = header[3].parse().unwrap();
        (num_nodes, source, sink)
    };

    let mut edges = Vec::new();
    for line in data[1..].iter() {
        let line: Vec<&str> = line.trim().split(' ').collect();
        let (from, to, capacity, gain) = (
            line[0].parse().unwrap(),
            line[1].parse().unwrap(),
            line[2].parse().unwrap(),
            line[3].parse().unwrap(),
        );
        edges.push((from, to, capacity, gain));
    }

    GraphInstancve {
        num_nodes,
        edges,
        source,
        sink,
    }
}
