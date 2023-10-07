use std::collections::VecDeque;
use std::fmt::Debug;

pub type Flow = i64;
pub const FLOW_MAX: Flow = Flow::MAX;

#[derive(Debug, Clone)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub flow: Flow,
    pub capacity: Flow,
}

#[derive(Debug)]
pub struct InsideEdge {
    pub to: usize,
    pub flow: Flow,
    pub capacity: Flow,
    pub rev: usize,
}

impl InsideEdge {
    #[inline]
    pub fn residual_capacity(&self) -> Flow {
        assert!(self.capacity >= self.flow);
        self.capacity - self.flow
    }
}

// CSR format
#[derive(Default)]
pub struct Graph {
    pub num_nodes: usize,
    pub num_edges: usize,
    pub edge_list: Vec<Edge>,

    pub start: Vec<usize>,
    pub inside_edge_list: Vec<InsideEdge>,

    pub excesses: Vec<Flow>,
    pub distance: Vec<usize>,
}

#[allow(dead_code)]
impl<'a> Graph {
    pub fn new() -> Self {
        Graph::default()
    }

    pub fn add_directed_edge(&mut self, from: usize, to: usize, capacity: Flow) -> Option<usize> {
        if capacity <= 0 as Flow {
            return None;
        }

        self.edge_list.push(Edge {
            from,
            to,
            flow: 0 as Flow,
            capacity,
        });
        self.num_nodes = self.num_nodes.max(from.max(to) + 1);
        self.num_edges += 1;
        Some(self.num_edges - 1)
    }

    pub fn get_directed_edge(&self, edge_index: usize) -> &Edge {
        &self.edge_list[edge_index]
    }

    pub fn build(&mut self) {
        let mut degree = vec![0; self.num_nodes];
        let mut edge_index = vec![usize::MAX; self.num_edges];
        let mut reverse_edge_index = vec![usize::MAX; self.num_edges];

        let mut tmp_inside_edge_list = Vec::with_capacity(2 * self.num_edges);
        for (i, e) in self.edge_list.iter().enumerate() {
            edge_index[i] = degree[e.from];
            degree[e.from] += 1;
            reverse_edge_index[i] = degree[e.to];
            degree[e.to] += 1;

            // from -> to
            tmp_inside_edge_list.push((
                e.from,
                InsideEdge {
                    to: e.to,
                    flow: 0 as Flow,
                    capacity: e.capacity,
                    rev: usize::MAX,
                },
            ));

            // to -> from
            tmp_inside_edge_list.push((
                e.to,
                InsideEdge {
                    to: e.from,
                    flow: e.capacity,
                    capacity: e.capacity,
                    rev: usize::MAX,
                },
            ));
        }

        // make graph
        self.excesses = vec![0 as Flow; self.num_nodes];
        self.distance = vec![0; self.num_nodes];
        self.start = vec![0; self.num_nodes + 1];
        self.inside_edge_list = (0..2 * self.num_edges)
            .map(|_| InsideEdge {
                to: 0,
                flow: 0 as Flow,
                capacity: 0 as Flow,
                rev: 0,
            })
            .collect();

        for (u, _) in tmp_inside_edge_list.iter() {
            self.start[u + 1] += 1;
        }
        for i in 1..=self.num_nodes {
            self.start[i] += self.start[i - 1];
        }

        let mut counter = self.start.clone();
        for (u, e) in tmp_inside_edge_list {
            self.inside_edge_list[counter[u]] = e;
            counter[u] += 1;
        }

        for (i, e) in self.edge_list.iter().enumerate() {
            edge_index[i] += self.start[e.from];
            reverse_edge_index[i] += self.start[e.to];
            self.inside_edge_list[edge_index[i]].rev = reverse_edge_index[i];
            self.inside_edge_list[reverse_edge_index[i]].rev = edge_index[i];
        }
    }

    pub fn neighbors(&'a self, u: usize) -> std::slice::Iter<'a, InsideEdge> {
        self.inside_edge_list[self.start[u]..self.start[u + 1]].iter()
    }

    pub fn push_flow(&mut self, u: usize, edge_index: usize, flow: Flow) {
        if flow == 0 as Flow {
            return;
        }
        let to = self.inside_edge_list[edge_index].to;
        let rev = self.inside_edge_list[edge_index].rev;

        // update flow
        self.inside_edge_list[edge_index].flow += flow;
        self.inside_edge_list[rev].flow -= flow;

        // update excess
        self.excesses[u] -= flow;
        self.excesses[to] += flow;
        assert!(
            self.inside_edge_list[edge_index].capacity >= self.inside_edge_list[edge_index].flow
                && self.inside_edge_list[edge_index].flow >= 0 as Flow
        );
        assert!(
            self.inside_edge_list[rev].capacity >= self.inside_edge_list[rev].flow
                && self.inside_edge_list[rev].flow >= 0 as Flow
        );
    }

    // O(n + m)
    // calculate distance from u to sink in residual network
    pub fn calculate_distance_to_sink(&mut self, sink: usize) -> Vec<usize> {
        let mut que = VecDeque::new();
        que.push_back(sink);
        let mut distance = vec![self.num_nodes; self.num_nodes];
        distance[sink] = 0;

        while let Some(u) = que.pop_front() {
            for edge in self.neighbors(u) {
                if edge.flow > 0 as Flow && distance[edge.to] > distance[u] + 1 {
                    distance[edge.to] = distance[u] + 1;
                    que.push_back(edge.to);
                }
            }
        }

        distance
    }

    #[inline]
    pub fn is_admissible_edge(&self, from: usize, to: usize) -> bool {
        self.distance[from] == self.distance[to] + 1
    }
}
