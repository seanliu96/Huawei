#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include <vector>

using namespace std;

const int MAX_V = 1010;
const int MAX_P_SIZE = 30;
const int inf = 0x3f3f3f3f;
const long long infll = 0x3f3f3f3f3f3f3f3f;
const int kmean_times = 5;
const double last_second = (45 >> 1) - 1, pm = 0.1, pc = 0.5, c1 = 1.0, c2 = 1.0, w = 0.9;


struct EdgeInfo {
    int v, w, c;
};

struct CustomerNodeInfo {
    int u, v, w;
};

struct Edge {
    int t, u, c, U;
    Edge *next, *pair;
};

class Particle {
public:
    Particle(int length=0);
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

class LiuXin {
public:
    void readtopo(char * topo[MAX_EDGE_NUM], int line_num);
    void spfa();
    vector<int> kmeans(int k);
    vector<vector<EdgeInfo> > graph;
    vector<CustomerNodeInfo> customer_nodes;
    vector<vector<int> > d;
    bool vis[MAX_V];
    int need_flow, node_num, edge_num, customer_num, server_cost;
};

class JinTao {
public:
    JinTao(LiuXin & lx);
    LiuXin *liuxin;
    void recover();
    void add_edge(int u, int v, int w, int c);
    void add_server(vector<int> &Q);
    long long costflow();
    void print_flow(vector<vector<int> > &node, vector<int> &flow);
private:
    int aug(int u, int m);
    bool modlabel();
    Edge epool[MAX_EDGE_NUM * 4 + 3000], *e[MAX_V];
    int psz, s, t;
    int flow, dist, d[MAX_V];
    long long cost;
    bool vis[MAX_V];
};

class HGAPSO {
public:
    HGAPSO(JinTao & jt, double pm, double pc, double c1, double c2, double w);
    void initial();
    vector<int> get_best();
    void addone(vector<int> & v);
    int run();
    int first_run();
private:
    Particle encode(vector<int> & vini);
    vector<int> decode(vector<double> & v);
    void GA_cross(Particle & s1, Particle & s2);
    void GA_mutation(Particle & s);
    void PSO_update(Particle & s);
    vector<Particle> p;
    Particle gbest;
    int l, unchanged_times;
    double GA_pm, GA_pc, PSO_c1, PSO_c2, PSO_w;
    JinTao *jintao;
};

template <class T>
void knuth_shuffle(vector<T> & v);

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename);

#endif
