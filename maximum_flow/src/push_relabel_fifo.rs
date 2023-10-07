use crate::graph::{Flow, Graph};
use std::collections::VecDeque;

#[derive(Default)]
pub struct PushRelabelFIFO {
    graph: Graph,

    active_nodes: VecDeque<usize>,
    current_edge: Vec<usize>,

    num_distance: Vec<usize>,
}

impl PushRelabelFIFO {
    pub fn new() -> Self {
        PushRelabelFIFO::default()
    }

    pub fn add_directed_edge(&mut self, from: usize, to: usize, capacity: Flow) -> Option<usize> {
        self.graph.add_directed_edge(from, to, capacity)
    }

    pub fn solve(&mut self, source: usize, sink: usize) -> Flow {
        self.graph.build();

        if source == sink || self.graph.num_nodes == 0 || self.graph.num_edges == 0 {
            return 0;
        }
        assert!(source < self.graph.num_nodes && sink < self.graph.num_nodes);

        self.pre_process(source, sink);

        while let Some(u) = self.active_nodes.pop_front() {
            if u == sink || self.graph.distance[u] >= self.graph.num_nodes {
                assert_ne!(u, source);
                continue;
            }
            self.discharge(u);
        }

        self.graph.excesses[sink]
    }

    fn pre_process(&mut self, source: usize, sink: usize) {
        self.current_edge = vec![0; self.graph.num_nodes];
        self.num_distance = vec![0; self.graph.num_nodes + 1];

        self.global_relabeling(sink);
        self.graph.distance[source] = self.graph.num_nodes;

        for u in 0..self.graph.num_nodes {
            self.num_distance[self.graph.distance[u]] += 1;
            self.current_edge[u] = self.graph.start[u];
        }

        for i in self.graph.start[source]..self.graph.start[source + 1] {
            let delta = self.graph.inside_edge_list[i].residual_capacity();
            self.graph.push_flow(source, i, delta);
        }

        for u in 0..self.graph.num_nodes {
            if u != source && u != sink && self.graph.excesses[u] > 0 {
                self.active_nodes.push_back(u);
            }
        }
    }

    fn discharge(&mut self, u: usize) {
        // push
        for i in self.current_edge[u]..self.graph.start[u + 1] {
            self.current_edge[u] = i;
            if self.graph.excesses[u] > 0 {
                self.push(u, i);
            }

            if self.graph.excesses[u] == 0 {
                return;
            }
        }
        self.current_edge[u] = self.graph.start[u];

        // relabel
        if self.num_distance[self.graph.distance[u]] == 1 {
            self.gap_relabeling(self.graph.distance[u]);
        } else {
            self.relabel(u);
        }

        if self.graph.excesses[u] > 0 {
            self.active_nodes.push_back(u);
        }
    }

    // push from u
    fn push(&mut self, u: usize, i: usize) {
        let to = self.graph.inside_edge_list[i].to;
        let delta = self.graph.excesses[u].min(self.graph.inside_edge_list[i].residual_capacity());
        if self.graph.is_admissible_edge(u, to) && delta > 0 {
            self.graph.push_flow(u, i, delta);
            if self.graph.excesses[to] == delta {
                self.active_nodes.push_back(to);
            }
        }
    }

    fn relabel(&mut self, u: usize) {
        self.num_distance[self.graph.distance[u]] -= 1;

        self.graph.distance[u] = self
            .graph
            .neighbors(u)
            .filter(|edge| edge.residual_capacity() > 0)
            .map(|edge| self.graph.distance[edge.to] + 1)
            .min()
            .unwrap()
            .min(self.graph.num_nodes);

        self.num_distance[self.graph.distance[u]] += 1;
    }

    // global relabeling heuristic
    // O(n + m)
    fn global_relabeling(&mut self, sink: usize) {
        self.graph.distance = self.graph.calculate_distance_to_sink(sink);
    }

    // gap relabeling heuristic
    fn gap_relabeling(&mut self, k: usize) {
        for u in 0..self.graph.num_nodes {
            if self.graph.distance[u] >= k {
                self.num_distance[self.graph.distance[u]] -= 1;
                self.graph.distance[u] = self.graph.distance[u].max(self.graph.num_nodes);
                self.num_distance[self.graph.distance[u]] += 1;
            }
        }
    }
}

#[cfg(test)]
mod test {
    use crate::graph::Flow;
    use crate::push_relabel_fifo::PushRelabelFIFO;
    use crate::test_utility::{
        read_expected, read_instance_aoj, read_instance_libreoj, GraphInstance,
    };
    use rstest::*;
    use std::path::PathBuf;

    #[rstest]
    fn aoj_grl_6_a(#[files("test_cases/AOJ_GRL_6_A/*.in")] path: PathBuf) {
        let mut expected_file_path = path.clone();
        expected_file_path.set_extension("out");

        let actual = execute(read_instance_aoj(&path));
        let expected = read_expected(&expected_file_path);
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn libreoj_101(#[files("test_cases/LibreOJ_101/*.in")] path: PathBuf) {
        let mut expected_file_path = path.clone();
        expected_file_path.set_extension("out");

        let actual = execute(read_instance_libreoj(&path));
        let expected = read_expected(&expected_file_path);
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn libreoj_127(#[files("test_cases/LibreOJ_127/*.in")] path: PathBuf) {
        let mut expected_file_path = path.clone();
        expected_file_path.set_extension("out");

        let actual = execute(read_instance_libreoj(&path));
        let expected = read_expected(&expected_file_path);
        assert_eq!(actual, expected);
    }

    fn execute(instance: GraphInstance) -> Flow {
        let mut solver = PushRelabelFIFO::new();
        for (from, to, capacity) in instance.edges {
            solver.add_directed_edge(from, to, capacity);
        }
        solver.solve(instance.source, instance.sink)
    }
}
