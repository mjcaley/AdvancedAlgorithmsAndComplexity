#include <algorithm>
#include <iostream>
#include <utility>
#include <queue>
#include <vector>

using std::vector;

/* This class implements a bit unusual scheme for storing edges of the graph,
 * in order to retrieve the backward edge for a given edge quickly. */
class FlowGraph {
public:
    struct Edge {
        int from, to, capacity, flow;
    };
    
private:
    /* List of all - forward and backward - edges */
    vector<Edge> edges;
    
    /* These adjacency lists store only indices of edges in the edges list */
    vector<vector<size_t> > graph;
    
public:
    explicit FlowGraph(size_t n): graph(n) {}
    
    void add_edge(int from, int to, int capacity) {
        /* Note that we first append a forward edge and then a backward edge,
         * so all forward edges are stored at even indices (starting from 0),
         * whereas backward edges are stored at odd indices in the list edges */
        Edge forward_edge = {from, to, capacity, 0};
        Edge backward_edge = {to, from, 0, 0};
        graph[from].push_back(edges.size());
        edges.push_back(forward_edge);
        graph[to].push_back(edges.size());
        edges.push_back(backward_edge);
    }
    
    size_t size() const {
        return graph.size();
    }
    
    const vector<size_t>& get_ids(int from) const {
        return graph[from];
    }
    
    const Edge& get_edge(size_t id) const {
        return edges[id];
    }
    
    void add_flow(size_t id, int flow) {
        /* To get a backward edge for a true forward edge (i.e id is even), we should get id + 1
         * due to the described above scheme. On the other hand, when we have to get a "backward"
         * edge for a backward edge (i.e. get a forward edge for backward - id is odd), id - 1
         * should be taken.
         *
         * It turns out that id ^ 1 works for both cases. Think this through! */
        edges[id].flow += flow;
        edges[id ^ 1].flow -= flow;
    }
};

vector<std::pair<int, int>> find_path(const FlowGraph& graph)
{
    auto source { 0 };
    auto sink { graph.size() -1 };
    
    using Node_Edge = std::pair<int, int>;
    vector<Node_Edge> parents(graph.size(), { -1, -1 });
    vector<bool> visited(graph.size());
    
    std::queue<int> Q;
    Q.push(source);
    
    while (!Q.empty())
    {
        auto current { Q.front() };
        Q.pop();
        visited[current] = true;
        
        // get the edge ids for current node
        auto edges { graph.get_ids(current) };
        for (auto edge_id : edges)
        {
            auto& edge { graph.get_edge(edge_id) };
            
            // break if we found the sink
            if (edge.to == sink && edge.capacity - edge.flow > 0)
            {
                parents[sink] = { current, edge_id };
                Q = std::queue<int>();  // empty queue
                break;
            }
            
            if (!visited[edge.to] && edge.capacity - edge.flow > 0)
            {
                Q.push(edge.to);
                parents[edge.to] = { current, edge_id };
            }
        }
    }
    
    
    vector<Node_Edge> path;
    
    auto parent { sink };
    while (parents[parent].second != -1)
    {
        path.emplace_back(parents[parent]);
        parent = parents[parent].first;
    }
    std::reverse(path.begin(), path.end());
    
    return path;
}

int min_flow(const FlowGraph& graph, const vector<std::pair<int, int>>& path)
{
    auto min = std::min_element(path.begin(),
                                path.end(),
                                [graph](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs)
                                {
                                    auto& lhs_edge = graph.get_edge(lhs.second);
                                    auto& rhs_edge = graph.get_edge(rhs.second);
                                    return (lhs_edge.capacity - lhs_edge.flow)
                                         < (rhs_edge.capacity - rhs_edge.flow);
                                });
    
    return graph.get_edge((*min).second).capacity - graph.get_edge((*min).second).flow;
}

void change_flow(FlowGraph& graph, const vector<std::pair<int, int>>& path, const int flow)
{
    for (auto& node_edge : path)
    {
        graph.add_flow(node_edge.second, flow);
    }
}

FlowGraph read_data() {
    int vertex_count, edge_count;
    std::cin >> vertex_count >> edge_count;
    FlowGraph graph(vertex_count);
    for (int i = 0; i < edge_count; ++i) {
        int u, v, capacity;
        std::cin >> u >> v >> capacity;
        graph.add_edge(u - 1, v - 1, capacity);
    }
    return graph;
}

int max_flow(FlowGraph& graph, int from, int to) {
    int flow = 0;
    /* your code goes here */
    
    while (true)
    {
        auto path = find_path(graph);
        if (path.empty())
        {
            break;
        }
        
        auto minimum_flow = min_flow(graph, path);
        change_flow(graph, path, minimum_flow);
        
        flow += minimum_flow;
    }
    
    return flow;
}

int main() {
    std::ios_base::sync_with_stdio(false);
    FlowGraph graph = read_data();
    
    std::cout << max_flow(graph, 0, graph.size() - 1) << "\n";
    return 0;
}
