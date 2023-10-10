use crate::graph::{Dist, Flow, InsideEdge, ScalingGraph, DIST_MAX, EPS, FLOW_MAX};

pub struct HighestGainPathMethod {
    graph: ScalingGraph,
}

#[allow(dead_code)]
impl HighestGainPathMethod {
    pub fn new(num_nodes: usize, epsilon: f64) -> Self {
        assert!(num_nodes > 0);
        assert!(epsilon > 0.0);

        HighestGainPathMethod {
            graph: ScalingGraph::new(num_nodes, epsilon),
        }
    }

    pub fn add_directed_edge(&mut self, from: usize, to: usize, capacity: Flow, gain: Flow) {
        self.graph.add_directed_edge(from, to, capacity, gain);
    }

    pub fn solve(&mut self, source: usize, sink: usize) -> Flow {
        self.graph.build();
        if self.graph.num_nodes == 0 || self.graph.num_edges == 0 {
            return 0 as Flow;
        }

        if !self.graph.is_lossy {
            match self.graph.calculate_distance_to_sink_with_negative_edge(sink) {
                Some(_distance_to_sink) => {}
                None => {
                    eprintln!("Error. flow generating cycle detected");
                    return 0.0;
                }
            }
        }

        self.graph.excesses[source] = FLOW_MAX;

        while self.graph.excesses[source] > EPS {
            if !self.argument_flow(source, sink) {
                break;
            }
        }

        self.graph.excesses[sink]
    }

    fn argument_flow(&mut self, source: usize, sink: usize) -> bool {
        match self.graph.find_shortest_path(source, sink) {
            None => false,
            Some(prev) => {
                // calculate delta and canonical labels
                let mut delta = FLOW_MAX;
                let mut canonical_labels = vec![FLOW_MAX; self.graph.num_nodes];
                canonical_labels[sink] = 1.0;

                let mut dist_to_sink = 0;
                let mut v = sink;
                while v != source {
                    // u -> v
                    let (u, i) = prev[v];
                    let edge = &self.graph.inside_edge_list[i];

                    dist_to_sink += edge.dist;

                    let label = self.graph.base.powf(dist_to_sink as f64);
                    canonical_labels[u] = label;

                    delta = delta.min(self.labeled_residual_capacity(u, edge, &canonical_labels));
                    v = u;
                }

                delta = delta.min(self.graph.excesses[source] / canonical_labels[source]);

                // update flow
                let mut v = sink;
                while v != source {
                    // u -> v
                    let (u, i) = prev[v];
                    self.graph.push_flow(u, i, delta, &canonical_labels);
                    v = u;
                }

                self.graph.excesses[source] -= canonical_labels[source] * delta;
                self.graph.excesses[sink] += delta;

                true
            }
        }
    }

    fn calculate_canonical_labels(&mut self, distance_to_sink: &Vec<Dist>, sink: usize) -> Vec<Flow> {
        let mut canonical_labels: Vec<Flow> = distance_to_sink
            .iter()
            .map(|&d| if d != DIST_MAX { self.graph.base.powf(d as Flow) } else { FLOW_MAX })
            .collect();
        canonical_labels[sink] = 1.0;
        canonical_labels
    }

    #[inline]
    fn labeled_residual_capacity(&self, u: usize, edge: &InsideEdge, labels: &Vec<f64>) -> Flow {
        edge.residual_capacity() / labels[u]
    }
}

#[cfg(test)]
mod tests {
    use super::HighestGainPathMethod;
    use super::{Flow, EPS};
    use crate::test_utilities::{read_expected, read_graph_instance};
    use rstest::*;
    use std::path::PathBuf;

    #[test]
    fn sample() {
        let epsilon = 0.01;
        let num_nodes = 8;
        let mut solver = HighestGainPathMethod::new(num_nodes, epsilon);
        solver.add_directed_edge(0, 1, 12.0, 0.7);
        solver.add_directed_edge(0, 2, 3.0, 0.9);
        solver.add_directed_edge(0, 3, 4.0, 0.8);

        solver.add_directed_edge(1, 4, 3.0, 0.5);
        solver.add_directed_edge(1, 5, 5.0, 0.8);

        solver.add_directed_edge(2, 1, 2.7, 1.0);
        solver.add_directed_edge(2, 3, 20.0 / 9.0, 0.9);
        solver.add_directed_edge(2, 5, 5.0, 0.7);

        solver.add_directed_edge(3, 5, 1.0, 1.0);
        solver.add_directed_edge(3, 6, 2.0, 0.7);

        solver.add_directed_edge(4, 7, 2.0, 0.5);

        solver.add_directed_edge(5, 4, 1.0, 0.5);
        solver.add_directed_edge(5, 6, 6.0, 0.7);
        solver.add_directed_edge(5, 7, 1.3, 1.0);

        solver.add_directed_edge(6, 7, 7.0, 1.0);

        let actual = solver.solve(0, 7);
        let expected = 7.363;
        assert!(expected * (1.0 - epsilon) <= actual && actual <= expected);
    }

    #[rstest]
    fn aoj_grl_6_a(#[files("test_cases/gain_random/*.in")] path: PathBuf) {
        let mut expected_file_path = path.clone();
        expected_file_path.set_extension("out");

        let epsilon: Flow = 0.01;
        let actual = get_result(&path, epsilon) as Flow;
        let expected = read_expected(&expected_file_path) as Flow;

        if expected == 0.0 {
            assert!(actual < EPS);
        } else {
            assert!(expected * (1.0 - epsilon) <= actual && actual <= expected);
        }
    }

    fn get_result(file_path: &PathBuf, epsilon: Flow) -> Flow {
        let instance = read_graph_instance(file_path);
        let mut solver = HighestGainPathMethod::new(instance.num_nodes, epsilon);

        for (from, to, capacity, gain) in instance.edges {
            solver.add_directed_edge(from, to, capacity, gain);
        }

        solver.solve(instance.source, instance.sink)
    }
}
