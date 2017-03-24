#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include <vector>

using namespace std;

const int MAX_V = 2010;
const int inf = 0x3f3f3f3f;
const long long infll = 0x3f3f3f3f3f3f3f3f;
//const double pm = 0.1, pc = 0.6, c1 = 1.0, c2 = 1.0, w = 0.9;
const double alpha = 0.995;
const double pm = 0.2, pc = 1.2, c1 = 2.0, c2 = 2.0, w = 0.9;

struct EdgeInfo {
    int v, c;
    EdgeInfo(int _v, int _c): v(_v), c(_c) {}
};

struct CustomerNodeInfo {
    int v, w;
    CustomerNodeInfo(int _v, int _w): v(_v), w(_w) {}
};

struct Edge {
    int t, u, c, U, C;
    Edge *next, *pair;
};

class Particle {
public:
    Particle(int length=0);
    Particle(int length, vector<int> & vi);
    vector<double> v;
    vector<double> v_best;
    vector<double> vp;
    long long cost_best;
    long long cost;
    friend bool operator== (const Particle & p1, const Particle & p2);
    friend bool operator< (const Particle & p1, const Particle & p2);
    friend bool operator<= (const Particle & p1, const Particle & p2);
    friend bool operator> (const Particle & p1, const Particle & p2);
    friend bool operator>= (const Particle & p1, const Particle & p2);
};

class Fuck {
public:
    void add_server(vector<int> & v);
    void add_edge(int u, int v, int w, int c);
    long long costflow();
    void print_flow(vector<vector<int> > &node, vector<int> &flow);
    void readtopo(char * topo[MAX_EDGE_NUM], int line_num);
    void spfa();
    void kmeans(int k, vector<int> & clusters);
    int need_flow, node_num, edge_num, customer_num;
    long long server_cost;
private:
    int aug(int u, int m);
    bool modlabel();
    vector<vector<EdgeInfo> > graph;
    vector<CustomerNodeInfo> customer_nodes;
    vector<vector<int> > d;
    bool vis[MAX_V];
    Edge epool[MAX_EDGE_NUM * 5], *e[MAX_V];
    int psz, s, t, psz_tmp;
    int flow, dist, D[MAX_V];
    long long cost;
};

class HGAPSO {
public:
    HGAPSO(Fuck & fuck);
    void get_best(vector<int> & server);
    void addone(vector<int> & v);
    int run();
    double initial(int max_p_size, int & best_size);
private:
    void decode(vector<double> & vd, vector<int> & vi);
    void GA_cross(Particle & s1, Particle & s2);
    void GA_mutation(Particle & s);
    void OBMA(Particle & s);
    void PSO_update(Particle & s);
    vector<Particle> p;
    vector<int> H;
    Particle gbest;
    int l, unchanged_times, iter;
    double GA_pm, GA_pc, PSO_c1, PSO_c2, PSO_w;
    Fuck *fuck;
};

template <class T>
void knuth_shuffle(vector<T> & v);

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename);

#endif
