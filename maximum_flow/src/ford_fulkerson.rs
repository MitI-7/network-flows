use crate::graph::{Flow, Graph};

#[derive(Default)]
pub struct FordFulkerson {
    graph: Graph,
}

impl FordFulkerson {
    pub fn new() -> Self {
        FordFulkerson::default()
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
            let mut used = vec![false; self.graph.num_nodes];
            let delta = self.dfs(source, sink, Flow::MAX, &mut used);
            if delta == 0 {
                return flow;
            }
            flow += delta;
        }
    }

    fn dfs(&mut self, u: usize, sink: usize, flow: Flow, visited: &mut Vec<bool>) -> Flow {
        if u == sink {
            return flow;
        }
        visited[u] = true;

        for i in self.graph.start[u]..self.graph.start[u + 1] {
            let to = self.graph.inside_edge_list[i].to;
            let residual_capacity = self.graph.inside_edge_list[i].residual_capacity();
            if visited[to] || residual_capacity == 0 {
                continue;
            }

            let delta = self.dfs(to, sink, flow.min(residual_capacity), visited);
            if delta > 0 {
                self.graph.push_flow(u, i, delta);
                return delta;
            }
        }
        0
    }
}

#[cfg(test)]
mod test {
    use crate::ford_fulkerson::{Flow, FordFulkerson};
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
        let mut solver = FordFulkerson::new();
        for (from, to, capacity) in instance.edges {
            solver.add_directed_edge(from, to, capacity);
        }
        solver.solve(instance.source, instance.sink)
    }
}
