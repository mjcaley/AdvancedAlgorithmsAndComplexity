#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <queue>

using std::vector;
using std::cin;
using std::cout;



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
    
    std::pair<vector<Edge>::const_iterator, vector<Edge>::const_iterator> get_edge_iterator() const
    {
        return std::make_pair(edges.begin(), edges.end());
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

FlowGraph build_graph(const vector<vector<bool>>& adj_matrix)
{
    int flights { static_cast<int>(adj_matrix.size()) };
    int crews { static_cast<int>(adj_matrix[0].size()) };
    
    FlowGraph graph(crews + flights + 2);
    
    int source { crews + flights };
    int sink { crews + flights + 1};
    
    // create edges from source to left-hand nodes
    for (int i { 0 }; i < flights; ++i)
    {
        graph.add_edge(source, i, 1);
    }
    // create edges from sink to right-hand nodes
    for (int i { flights }; i < flights + crews; ++i)
    {
        graph.add_edge(i, sink, 1);
    }
    
    // process adjacency matrix to create middle edges
    for (int i { 0 }; i < flights; ++i)
    {
        for (int j { 0 }; j < crews; ++j)
        {
            if (adj_matrix[i][j] == 1)
            {
                graph.add_edge(i, flights + j, 1);
            }
        }
    }
    
    return graph;
}

vector<std::pair<int, int>> find_path(const FlowGraph& graph, int from, int to)
{
    auto source { from };
    auto sink { to };
    
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

int max_flow(FlowGraph& graph, int from, int to) {
    int flow = 0;
    
    while (true)
    {
        auto path = find_path(graph, from, to);
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

vector<int> find_matching(const FlowGraph& graph, int num_flight, int num_crew)
{
    vector<int> matching(num_flight, -1);
    
    auto iters { graph.get_edge_iterator() };
    for (; iters.first != iters.second; ++iters.first)
    {
        if (iters.first->flow > 0
            && iters.first->from != graph.size() - 2
            && iters.first->to != graph.size() - 1)
        {
            matching[iters.first->from] = iters.first->to - num_flight;
        }
    }
    
    return matching;
}


class MaxMatching {
public:
    void Solve() {
        vector<vector<bool>> adj_matrix = ReadData();
        vector<int> matching = FindMatching(adj_matrix);
        WriteResponse(matching);
    }
    
private:
    vector<vector<bool>> ReadData() {
        int num_left, num_right;
        cin >> num_left >> num_right;
        vector<vector<bool>> adj_matrix(num_left, vector<bool>(num_right));
        for (int i = 0; i < num_left; ++i)
            for (int j = 0; j < num_right; ++j) {
                int bit;
                cin >> bit;
                adj_matrix[i][j] = (bit == 1);
            }
        return adj_matrix;
    }
    
    void WriteResponse(const vector<int>& matching) {
        for (int i = 0; i < matching.size(); ++i) {
            if (i > 0)
                cout << " ";
            if (matching[i] == -1)
                cout << "-1";
            else
                cout << (matching[i] + 1);
        }
        cout << "\n";
    }
    
    vector<int> FindMatching(const vector<vector<bool>>& adj_matrix) {
        auto graph = build_graph(adj_matrix);
        max_flow(graph, graph.size() - 2, graph.size() - 1);
        auto matching { find_matching(graph, adj_matrix.size(), adj_matrix[0].size()) };
        
        return matching;
    }
};

int main() {
    std::ios_base::sync_with_stdio(false);
    MaxMatching max_matching;
    max_matching.Solve();
    return 0;
}
