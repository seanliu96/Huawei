#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
using namespace std;

struct EdgeInfo {
    int v, w, c;
};

struct CustomerNodeInfo {
    int u, v, w;
};

struct ServerType {
    ServerType(){}
    ServerType(int _flow, int _cost) {
        cost = _cost;
        flow = _flow;
    }
    int cost, flow;
};

const int inf = 0x3f3f3f3f;
const long long infll = 0x3f3f3f3f3f3f3f3f;
long long cost = 0;
vector<int> servers_flow;
vector<vector<EdgeInfo> > graph;
vector<CustomerNodeInfo> customer_nodes;
vector<ServerType> Server;  
vector<int> node_cost;
int node_num, edge_num, customer_num, server_cost, need_flow;

bool readtopo(char * filename);
long long judge(char * filename);

int main(int argc, char * argv[]) {
    if (argc != 3) {
        cout << "error: give me the cast file and output file\n";
        return 1;
    }
    if(readtopo(argv[1]))
        cout << judge(argv[2]) << endl;
}

bool readtopo(char * filename) {
    int u, v, c, w;
    fstream fin(filename, ios::in);
    if (!fin.is_open()) {
        cout << "error: fail to open" << filename << endl;
        return false;
    }
    fin >> node_num >> edge_num >> customer_num;
    graph.resize(node_num, vector<EdgeInfo>());
    customer_nodes.resize(customer_num);
    node_cost.resize(node_num);
    servers_flow.resize(node_num, 0);
    fin >> u >> w >> c;
    Server.push_back(ServerType(w, c));
    while (fin >> u) {
        if (u == 0)
            break;
        fin >> w >> c;
        Server.push_back(ServerType(w, c));
    }
    fin >> c;
    node_cost[u] = c;
    for (int i = 1; i < node_num; ++i) {
        fin >> u >> c;
        node_cost[u] = c;
    }
    for (int i = 0; i < edge_num; ++i) {
        fin >> u >> v >> w >> c;
        graph[u].push_back((EdgeInfo){v, w, c});
        graph[v].push_back((EdgeInfo){u, w, c});
    }
    need_flow = 0;
    for (int i = 0; i < customer_num; ++i) {
        fin >> u >> v >> w;
        customer_nodes[i] = (CustomerNodeInfo){u, v, w};
        need_flow += w;
    }
    fin.close();
    return true;
}

long long judge(char * filename) {
    long long cost = 0;
    int size, n, flow;
    fstream fin(filename, ios::in);
    if (!fin.is_open()) {
        cost = infll;
        return cost;
    }
    vector<int> line;
    string buf;
    getline(fin, buf);
    stringstream ss(buf);
    ss >> size;
    getline(fin, buf);
    for (int i = 0; i < size; ++i) {
        line.clear();
        getline(fin, buf);
        stringstream sss(buf);
        while (sss >> n) {
            line.push_back(n);
        }
        flow = line[line.size() - 2];
        if (!servers_flow[line[0]]) {
            servers_flow[line[0]] = Server[line[line.size() - 1]].flow;
            cost += Server[line[line.size() - 1]].cost;
            cost += node_cost[line[0]];
        } 
        for (int j = 0; j < line.size() - 4; ++j) {
            bool find = false;
            for (int k = 0; k < graph[line[j]].size(); ++k) {
                if (graph[line[j]][k].v == line[j+1]) {
                    find = true;
                    graph[line[j]][k].w -= flow;
                    cost += graph[line[j]][k].c * flow;
                    if (graph[line[j]][k].w < 0) {
                        cost = infll;
                        fin.close();
                        return cost;
                    }
                    break;
                }
            }
            if (!find) {
                cost = infll;
                fin.close();
                return cost;
            }
        }
        customer_nodes[line[line.size() - 3]].w -= flow;
        servers_flow[line[0]] -= flow;

    }
    for (int i = 0; i < customer_nodes.size(); ++i) {
        if (customer_nodes[i].w > 0) {
            cost = infll;
            fin.close();
            return cost;
        }
    }

    for (int i = 0; i < node_num; ++i) {
        if (servers_flow[i] < 0) {
            cost = infll;
            fin.close();
            return cost;
        }
    }
    fin.close();
    return cost;
}
