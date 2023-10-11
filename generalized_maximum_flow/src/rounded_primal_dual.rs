use crate::graph::{Dist, Flow, InsideEdge, ScalingGraph, DIST_MAX, EPS, FLOW_MAX};
use std::collections::VecDeque;

pub struct RoundedPrimalDual {
    base: Flow,
    pub graph: ScalingGraph,
    pub excesses: Vec<Flow>,
    labels: Vec<Flow>,

    // maximum flow(dinic)
    iter: Vec<usize>,
    level: Vec<isize>,
}

#[allow(dead_code)]
impl RoundedPrimalDual {
    pub fn new(num_nodes: usize, epsilon: Flow) -> Self {
        let base = (1.0 + epsilon).powf(1.0 / num_nodes as Flow) as Flow;
        Self::new_with_base(base)
    }

    pub fn new_with_base(base: Flow) -> Self {
        assert!(base > 1.0);

        RoundedPrimalDual {
            base,
            graph: ScalingGraph::new_with_base(base),
            excesses: Vec::new(),
            labels: Vec::new(),

            iter: Vec::new(),
            level: Vec::new(),
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

        self.excesses = vec![0.0; self.graph.num_nodes];
        self.labels = vec![0.0; self.graph.num_nodes];

        self.excesses[source] = FLOW_MAX;

        if !self.graph.is_lossy {
            match self.graph.calculate_distance_to_sink_with_negative_edge(sink) {
                Some(distance_to_sink) => {
                    self.update_labels(&distance_to_sink, sink);
                }
                None => {
                    eprintln!("Error. flow generating cycle detected");
                    return 0.0;
                }
            }
        }

        while self.excesses[source] > EPS {
            if !self.argument_flow(source, sink) {
                break;
            }
        }

        self.excesses[sink]
    }

    fn argument_flow(&mut self, source: usize, sink: usize) -> bool {
        let distance_to_sink = self.graph.calculate_distance_to_sink(sink);
        self.update_labels(&distance_to_sink, sink);

        // no augmenting path from source
        if self.labels[source] == FLOW_MAX {
            return false;
        }

        // maximum flow
        while self.excesses[source] > EPS as Flow {
            self.bfs(source);
            if self.level[sink] < 0 {
                break;
            }
            self.iter = (0..self.graph.num_nodes).map(|u| self.graph.start[u]).collect();
            while self.excesses[source] > EPS {
                let flow = self.dfs(source, sink, self.excesses[source]);
                if flow <= 0 as Flow {
                    break;
                }
                self.excesses[source] -= flow * self.labels[source];
                self.excesses[sink] += flow;
            }
        }

        true
    }

    fn update_labels(&mut self, distance_to_sink: &Vec<Dist>, sink: usize) {
        self.labels = distance_to_sink
            .iter()
            .map(|&d| if d != DIST_MAX { self.base.powf(d as Flow) } else { FLOW_MAX })
            .collect();
        self.labels[sink] = 1.0;
    }

    fn bfs(&mut self, source: usize) {
        self.level = vec![-1; self.graph.num_nodes];
        let mut que = VecDeque::new();
        self.level[source] = 0;
        que.push_back(source);

        while let Some(u) = que.pop_front() {
            for i in self.graph.start[u]..self.graph.start[u + 1] {
                let edge = &self.graph.inside_edge_list[i];
                // for edge in self.graph.graph[u].iter() {
                if edge.residual_capacity() > 0.0 && self.level[edge.to] < 0 && self.reduced_cost(u, edge) == 0 {
                    self.level[edge.to] = self.level[u] + 1;
                    que.push_back(edge.to);
                }
            }
        }
    }

    fn dfs(&mut self, u: usize, sink: usize, flow: Flow) -> Flow {
        if u == sink {
            return flow;
        }
        for i in self.iter[u]..self.graph.start[u + 1] {
            self.iter[u] = i;
            let edge = &self.graph.inside_edge_list[i];
            if edge.residual_capacity() > 0.0 && self.level[u] < self.level[edge.to] && self.reduced_cost(u, &edge) == 0 {
                let d = flow.min(self.labeled_residual_capacity(u, edge));
                let f = self.dfs(edge.to, sink, d);
                if f > 0.0 {
                    self.graph.push_flow(u, i, f, &self.labels);
                    return f;
                }
            }
        }

        0.0
    }

    #[inline]
    fn labeled_residual_capacity(&self, u: usize, edge: &InsideEdge) -> Flow {
        edge.residual_capacity() / self.labels[u]
    }

    #[inline]
    fn reduced_cost(&self, u: usize, edge: &InsideEdge) -> Dist {
        edge.dist - self.graph.potentials[u] + self.graph.potentials[edge.to]
    }
}

#[cfg(test)]
mod tests {
    use super::{Flow, RoundedPrimalDual, EPS};
    use crate::test_utilities::{read_expected, read_graph_instance};
    use rstest::*;
    use std::path::PathBuf;

    #[test]
    fn sample() {
        let epsilon: Flow = 0.01;
        let num_nodes = 8;
        let mut solver = RoundedPrimalDual::new(num_nodes, epsilon);
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
        let expected = 7.363 as Flow;
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
        let mut solver = RoundedPrimalDual::new(instance.num_nodes, epsilon);

        for (from, to, capacity, gain) in instance.edges {
            solver.add_directed_edge(from, to, capacity, gain);
        }

        solver.solve(instance.source, instance.sink)
    }
}
