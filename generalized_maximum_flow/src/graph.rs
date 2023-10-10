use std::cmp::Reverse;
use std::collections::BinaryHeap;

pub type Flow = f64;
pub type Dist = i32;

pub const DIST_MAX: Dist = Dist::MAX / 2;
pub const FLOW_MAX: Flow = Flow::MAX / 2.0;
pub const EPS: Flow = Flow::EPSILON;

#[derive(Debug, Clone)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub flow: Flow,
    pub capacity: Flow,
    pub gain: Flow,
}

#[derive(Debug, Clone, Default)]
pub struct InsideEdge {
    pub to: usize,
    pub flow: Flow,
    pub capacity: Flow,
    pub dist: Dist,
    pub rev: usize,
}

impl InsideEdge {
    #[inline]
    pub fn residual_capacity(&self) -> Flow {
        // assert!(self.capacity >= self.flow);
        self.capacity - self.flow
    }
}

pub struct ScalingGraph {
    pub num_nodes: usize,
    pub num_edges: usize,
    pub base: Flow,
    edge_list: Vec<Edge>,
    pub is_lossy: bool,

    pub start: Vec<usize>,
    pub inside_edge_list: Vec<InsideEdge>,

    pub excesses: Vec<Flow>,
    pub potentials: Vec<Dist>,
}

#[allow(dead_code)]
impl<'a> ScalingGraph {
    pub fn new(num_nodes: usize, epsilon: Flow) -> Self {
        let base = (1.0 + epsilon).powf(1.0 / num_nodes as Flow) as Flow;
        Self::new_with_base(base)
    }

    pub fn new_with_base(base: Flow) -> Self {
        ScalingGraph {
            num_nodes: 0,
            num_edges: 0,
            base,
            edge_list: Vec::new(),
            is_lossy: true,

            start: Vec::new(),
            inside_edge_list: Vec::new(),

            excesses: Vec::new(),
            potentials: Vec::new(),
        }
    }

    pub fn add_directed_edge(&mut self, from: usize, to: usize, capacity: Flow, gain: Flow) -> Option<usize> {
        if gain <= 0.0 {
            eprintln!("warning. gain needs to be greater than 0");
            return None;
        }
        if capacity <= 0.0 {
            eprintln!("warning. capacity needs to be greater than 0");
            return None;
        }

        if gain > 1.0 {
            self.is_lossy = false;
        }

        self.edge_list.push(Edge {
            from,
            to,
            flow: 0 as Flow,
            capacity,
            gain,
        });
        self.num_nodes = self.num_nodes.max(from.max(to) + 1);
        self.num_edges += 1;

        Some(self.num_edges - 1)
    }

    pub fn get_directed_edge(&self, edge_index: usize) -> &Edge {
        &self.edge_list[edge_index]
    }

    pub fn neighbors(&'a self, u: usize) -> std::slice::Iter<'a, InsideEdge> {
        self.inside_edge_list[self.start[u]..self.start[u + 1]].iter()
    }

    #[inline]
    pub fn push_flow(&mut self, u: usize, i: usize, flow: Flow, labels: &Vec<Flow>) {
        let to = self.inside_edge_list[i].to;
        let rev = self.inside_edge_list[i].rev;
        self.inside_edge_list[i].flow += flow * labels[u];
        self.inside_edge_list[rev].flow -= flow * labels[to];

        if self.inside_edge_list[i].flow > self.inside_edge_list[i].capacity {
            self.inside_edge_list[i].flow = self.inside_edge_list[i].capacity;
            self.inside_edge_list[rev].flow = 0.0;
        }

        if self.inside_edge_list[rev].flow < 0.0 {
            self.inside_edge_list[rev].flow = 0.0;
            self.inside_edge_list[i].flow = self.inside_edge_list[i].capacity;
        }

        if self.inside_edge_list[i].residual_capacity() <= EPS || self.inside_edge_list[rev].flow <= EPS {
            self.inside_edge_list[i].flow = self.inside_edge_list[i].capacity;
            self.inside_edge_list[rev].flow = 0.0;
        }
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

            // gain scaling
            let c = e.gain.log(self.base).floor();
            let scaled_gain = self.base.powf(c);
            let dist = -c as Dist; // TODO: check over flow

            // from -> to
            tmp_inside_edge_list.push((
                e.from,
                InsideEdge {
                    to: e.to,
                    flow: 0 as Flow,
                    capacity: e.capacity,
                    dist,
                    rev: usize::MAX,
                },
            ));

            // to -> from
            tmp_inside_edge_list.push((
                e.to,
                InsideEdge {
                    to: e.from,
                    flow: e.capacity * scaled_gain,
                    capacity: e.capacity * scaled_gain,
                    dist: -dist,
                    rev: usize::MAX,
                },
            ));
        }

        // make graph
        self.start = vec![0; self.num_nodes + 1];
        self.excesses = vec![0 as Flow; self.num_nodes];
        self.potentials = vec![0; self.num_nodes];
        self.inside_edge_list = vec![Default::default(); 2 * self.num_edges];

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
}

impl ScalingGraph {
    // find shortest path from source to sink & update potentials
    pub fn find_shortest_path(&mut self, source: usize, sink: usize) -> Option<Vec<(usize, usize)>> {
        let mut prev = vec![(self.num_nodes, self.num_nodes); self.num_nodes];

        let mut visited = vec![false; self.num_nodes];
        let mut distance = vec![DIST_MAX; self.num_nodes];
        distance[source] = 0;

        let mut heap = BinaryHeap::new();
        heap.push((Reverse(0), source));
        while let Some((d, u)) = heap.pop() {
            if visited[u] {
                continue;
            }
            visited[u] = true;
            if u == sink {
                break;
            }

            let st = self.start[u];
            let en = self.start[u + 1];
            for (i, e) in self.inside_edge_list[st..en].iter().enumerate() {
                if e.residual_capacity() < EPS {
                    continue;
                }
                if visited[e.to] {
                    continue;
                }

                let dist = e.dist + self.potentials[u] - self.potentials[e.to];
                assert!(dist >= 0);
                let dist = dist;

                let new_dist = d.0 + dist;
                if new_dist < distance[e.to] {
                    distance[e.to] = new_dist;
                    prev[e.to] = (u, st + i);
                    heap.push((Reverse(new_dist), e.to));
                }
            }
        }

        // update potentials
        for u in 0..self.num_nodes {
            if visited[u] {
                self.potentials[u] += distance[u] as Dist - distance[sink] as Dist;
            }
        }

        if !visited[sink] {
            return None;
        }

        Some(prev)
    }
}
