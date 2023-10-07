use crate::graph::{Flow, Graph};
use std::collections::VecDeque;

#[derive(Default)]
pub struct Dinic {
    graph: Graph,
    current_edge: Vec<usize>,
    level: Vec<isize>,
}

impl Dinic {
    pub fn new() -> Self {
        Dinic::default()
    }

    pub fn add_directed_edge(&mut self, from: usize, to: usize, capacity: Flow) -> Option<usize> {
        self.graph.add_directed_edge(from, to, capacity)
    }

    pub fn solve(&mut self, source: usize, sink: usize) -> Flow {
        self.graph.build();
        if source == sink || self.graph.num_nodes == 0 || self.graph.num_edges == 0 {
            return 0;
        }

        let mut flow = 0;
        loop {
            self.bfs(source);
            if self.level[sink] < 0 {
                return flow;
            }

            self.current_edge = (0..self.graph.num_nodes)
                .map(|u| self.graph.start[u])
                .collect();
            loop {
                let delta = self.dfs(source, sink, Flow::MAX);
                if delta == 0 {
                    break;
                }
                flow += delta;
            }
        }
    }

    fn bfs(&mut self, source: usize) {
        self.level = vec![-1; self.graph.num_nodes];
        let mut que = VecDeque::new();
        self.level[source] = 0;
        que.push_back(source);

        while let Some(u) = que.pop_front() {
            for edge in self.graph.neighbors(u) {
                if edge.residual_capacity() > 0 && self.level[edge.to] < 0 {
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

        for i in self.current_edge[u]..self.graph.start[u + 1] {
            self.current_edge[u] = i;
            let edge = &self.graph.inside_edge_list[i];
            let to = edge.to;
            let residual_capacity = edge.residual_capacity();

            if residual_capacity > 0 && self.level[u] + 1 == self.level[to] {
                let d = self.dfs(to, sink, flow.min(residual_capacity));
                if d > 0 {
                    self.graph.push_flow(u, i, d);
                    return d;
                }
            }
        }
        self.current_edge[u] = self.graph.start[u + 1];

        0
    }
}

#[cfg(test)]
mod test {
    use crate::dinic::{Dinic, Flow};
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

    fn execute(instance: GraphInstance) -> Flow {
        let mut solver = Dinic::new();
        for (from, to, capacity) in instance.edges {
            solver.add_directed_edge(from, to, capacity);
        }
        solver.solve(instance.source, instance.sink)
    }
}
