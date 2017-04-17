#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include <vector>

using namespace std;

const int MAX_V = 2010;
const int inf = 0x3f3f3f3f;
const long long infll = 0x3f3f3f3f3f3f3f3f;
const double c1 = 1.0, c2 = 1.6, w = 0.9;

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

class Fuck {
public:
    void add_server(vector<int> & v);
    inline void add_edge(int u, int v, int w, int c);
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
    int vis[MAX_V];
    Edge epool[MAX_EDGE_NUM * 5], *e[MAX_V];
    int psz, s, t, psz_tmp;
    int flow, dist, D[MAX_V];
    long long cost;
};

class Particle {
public:
    Particle(int length=0);
    Particle(int length, vector<int> & vi, Fuck* & fuck);
    vector<double> v;
    vector<double> v_best;
    vector<double> vp;
    long long cost_best;
    long long cost;
};

class XJBS {
public:
    XJBS(Fuck & fuck);
    inline void get_best(vector<int> & server);
    inline void addone(vector<int> & v);
    void run1();
    void run2();
    void initial();
    inline void reproduction();
    int max_p_size;
    Particle gbest;
private:
    inline void decode(vector<double> & vd, vector<int> & vi);
    inline void GA_cross(Particle & s1, Particle & s2);
    inline void OBMA(Particle & s);
    inline void PSO_update(Particle & s);
    inline void updateone(Particle & s);
    vector<Particle> p;
    int l, cnt;
    double PSO_c1, PSO_c2, PSO_w;
    Fuck *fuck;
};

template <class T>
inline void knuth_shuffle(vector<T> & v);

bool cmp(const Particle & p1, const Particle & p2);

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename);

#endif
