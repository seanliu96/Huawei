#include "deploy.h"
#include <stdio.h>
#include <deque>
#include <queue>
#include <vector>
#include <limits.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstring>
#include <ctime>
#include <utility>
#include <cmath>

using namespace std;
//你要完成的功能总入口

clock_t run_second, last_second = (90 - 1.1) * CLOCKS_PER_SEC;

Fuck fuck;

vector<int> server, best_server;

vector<ServerType> Server;

int node_cost[MAX_V];
vector<vector<int> > node;
vector<int> flow;

FlowCost flow_cost[10005];
int node_server_id[MAX_V];


void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    fuck.readtopo(topo, line_num);
    fuck.spfa();
    XJBS xjbs(fuck);
    int best_index = fuck.customer_num;
    int kmean_times = 2;
    int up = fuck.customer_num * 0.8;
    int down = fuck.customer_num * 0.2;
    if (!up) up = 1;
    if (!down) down = 1;
    int block_size = (up - down) / 60 + 1;
    //fuck.kmeans(1, server);
    //xjbs.addone(server);
    fuck.kmeans(best_index, best_server);
    long long best_cost = fuck.costflow();
    
    for (int i = down; i <= up; i += block_size) {
        for (int j = 0; j < kmean_times ; ++j) {
            fuck.kmeans(i, server);
            fuck.add_server(server);
            long long cost = fuck.costflow();
            if (cost < best_cost) {
                best_cost = cost;
                best_index = i;
                best_server.swap(server);
            }
        }
    }
    xjbs.addone(best_server);
    xjbs.initial();
    run_second = last_second * 0.4;
    while (clock() < run_second) {
        xjbs.run1();
    }
    xjbs.reproduction();
    run_second += last_second * 0.6;
    while (clock() < run_second) {
        xjbs.run2();
    }
    xjbs.get_best(best_server, best_cost);

    fuck.add_server(best_server);
    fuck.costflow();
    fuck.print_flow(node, flow);
    int node_size = node.size();
    char * topo_file = new char[node_size * MAX_V * 5];
    topo_file[0] = '\0';
    char line[MAX_V * 5];
    char tmp[100];
    sprintf(line, "%d\n\n", node_size);
    strcat(topo_file, line);
    for (int i = 0; i < node_size; ++i) {
        line[0] = '\0';
        int node_size_1 = node[i].size() - 1;
        for (int j = 0; j < node_size_1; ++j) {
            sprintf(tmp, "%d ", node[i][j]);
            strcat(line, tmp);
        }
        sprintf(tmp, "%d ", node[i][node_size_1] - fuck.node_num);
        strcat(line, tmp);
        sprintf(tmp, "%d ", (int)flow[i]);
        strcat(line, tmp);
        sprintf(tmp, "%d\n", (int)node_server_id[node[i][0]]);
        strcat(line, tmp);
        strcat(topo_file, line);
    }
    write_result(topo_file, filename);
    delete []topo_file;

}

template <class T>
inline void knuth_shuffle(vector<T> & v) {
    int i = v.size() - 1, j = 0;
    while (i >= 0) {
        j = rand() % (i + 1);
        swap(v[i], v[j]);
        --i;
    }
}

void Fuck::readtopo(char * topo[MAX_EDGE_NUM], int line_num) { 
    int line = 0;
    int u, v, c, w;
    if (line < line_num)
        sscanf(topo[line], "%d %d %d", &node_num, &edge_num, &customer_num);
    s = node_num + customer_num;
    t = s + 1;
    psz = 0;
    memset(e, 0, sizeof(e));
    graph.resize(node_num, vector<EdgeInfo>());
    d.resize(node_num, vector<int>(node_num, inf));
    line += 2;
    Max_flow = 0;
    while (true) {
        if (strcmp("\r\n", topo[line]) == 0) break;
        sscanf(topo[line], "%d %d %d", &u, &w, &c);
        Server.push_back(ServerType(w, c));
        Max_flow = max(Max_flow, w);
        ++line;
    }
    for (int i = 0; i <= Max_flow; ++i) {
        flow_cost[i].cost = inf;
        for (int j = 0; j < Server.size(); ++j) 
        if (Server[j].flow >= i && Server[j].cost < flow_cost[i].cost) {
            flow_cost[i].cost = Server[j].cost;
            flow_cost[i].id = j;
        }
    }
    ++line;
    while (true) {
        if (strcmp("\r\n", topo[line]) == 0) break;
        sscanf(topo[line], "%d %d", &u, &c);
        node_cost[u] = c;
        ++line;
    }
    ++line;
    for (int i = 0; i < edge_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d %d", &u, &v, &w, &c);
        graph[u].emplace_back(v, c);
        graph[v].emplace_back(u, c);
        add_edge(u, v, w, c);
        add_edge(v, u, w, c);
    }
    ++line;
    need_flow = 0;
    for (int i = 0; i < customer_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d", &u, &v, &w);
        customer_nodes.emplace_back(v, w);
        add_edge(v, u + node_num, w, 0);
        add_edge(u + node_num, t, w, 0);
        need_flow += w;
    }
    psz_tmp = psz;
}

void Fuck::spfa() {
    for(int s = 0; s < node_num; ++s) {
        d[s][s] = 0;
        deque<int> q;
        memset(vis,0,sizeof(vis));
        vis[s] = 1;
        q.push_back(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            vis[u] = 0;
            for (unsigned int i = 0; i < graph[u].size(); ++i) {
                int v = graph[u][i].v;
                int dis =  d[s][u] + graph[u][i].c;
                if (dis < d[s][v])
                { 
                    d[s][v] = dis;
                    if (!vis[v]) {
                        vis[v] = 1;
                        if (q.size () && d[s][v] < d[s][q[0]])
                            q.push_front(v);
                        else
                            q.push_back(v);
                    }
                }
            }
        }
    }
}

void Fuck::kmeans(int k, vector<int> & clusters) {
    clusters.resize(k);
    memset(vis, -1, sizeof(vis));
    vector<vector<int> > kmean_node(k);
    int min_dist, min_index;
    knuth_shuffle(customer_nodes);
    for (int i = 0; i < k; ++i) {
        clusters[i] = customer_nodes[i].v;
    }
    while (1) {
        for (int i = 0; i < k; ++i)
            kmean_node[i].clear();
        bool update = 0;
        for (int i = 0; i < customer_num; ++i) {
            min_dist = inf;
            min_index = 0;
            for (int j = 0; j < k; ++j) {
                if(d[clusters[j]][customer_nodes[i].v] < min_dist) {
                    min_dist = d[clusters[j]][customer_nodes[i].v];
                    min_index = j;
                }
            }
            if (vis[i] != min_index) {
                update = 1;
                vis[i] = min_index;
            }
            kmean_node[vis[i]].push_back(i);
        }
        if (!update) 
            break;
        for (int j = 0; j < k; ++j) {
            min_dist = inf;
            min_index = 0;
            for (int l = 0; l < node_num; ++l) {
                int dist = 0;
                for (unsigned int i = 0; i < kmean_node[j].size(); ++i) {
                    dist += d[l][customer_nodes[kmean_node[j][i]].v];
                }
                if (dist < min_dist) {
                    min_index = l;
                    min_dist = dist;
                }
            }
            clusters[j] = min_index;
        }
    }
}

void Fuck::add_server(vector<int> & v) {
    if (psz != psz_tmp) {
        psz = psz_tmp;
        for (Edge *j = e[s]; j; j = j->next) {
            int x = j->t;
            e[x] = e[x]->next;
        }
        e[s] = 0;
        Edge *j = epool + psz;
        for (Edge *i = epool; i < j; ++i) {
            i->u = i->U;
            i->c = i->C;
        }
    }
    for (unsigned int i = 0; i < v.size(); ++i) {
        add_edge(s, v[i], Max_flow, 0);
    }
}

inline void Fuck::add_edge(int u, int v, int w, int c) {
    Edge *e1 = epool + psz++, *e2 = epool + psz++;
    *e1 = (Edge){v, w, c, w, c, e[u], e2}, e[u] = e1;
    *e2 = (Edge){u, 0, -c, 0, -c, e[v], e1}, e[v] = e2;
}

void Fuck::print_flow(vector<vector<int> > & node, vector<int> &flow) {
    node.clear();
    flow.clear();
    for (Edge *i = e[s]; i; i = i->next) {
        node_server_id[i->t] = flow_cost[i->U - i->u].id;
    }
    while (1) {
        vector<int> Tmp;
        int u = s;
        int S = inf;
        while (u != t) {
            bool flag=0;
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    S = min(S, i->U - i->u);
                    u = v;
                    flag = 1;
                    break;
                }
            }
            if (!flag) break;
        }
        if (u != t) break;
        u = s;
        flow.push_back(S);
        while (u != t) {
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    i->u += S;
                    u = v;
                    break;
                }
            }
            if (u != t) Tmp.push_back(u);
        }
        node.push_back(Tmp);
    }
}

int Fuck::aug(int u, int m) {
    if (u == t)
        return cost += (long long)dist * m, m;
    int d = m;
    vis[u] = 1;
    for (Edge *i = e[u]; i; i = i->next) {
        if (i->u && !i->c && !vis[i->t]) {
            int f = aug(i->t, min(d, i->u));
            i->u -= f;
            i->pair->u += f;
            d -= f;
            if (!d)
                return m;
        }
    }
    return m - d;
}

bool Fuck::modlabel() {
    deque <int> q;
    memset(vis , 0, sizeof(vis));
    memset(D, 0x3f , sizeof(D));
    q.push_back(s);
    D[s] = 0;
    vis[s] = 1;
    while (!q.empty ()) {
        int u = q.front ();
        q.pop_front ();
        for (Edge *i = e[u]; i; i = i->next) {
            int v = i->t;
            int dis = D[u] + i->c;
            if (i->u && dis < D[v]) {
                D[v] = dis;
                if (!vis[v]) {
                    vis[v] = 1;
                    if (q.size () && D[v] < D[q[0]])
                        q.push_front(v);
                    else
                        q.push_back(v);
                }
            }
        }
        vis[u] = 0;
    }
    for (Edge *i = epool; i < epool + psz; ++i) {
        i->c -= D[i->t] - D[i->pair->t];
    }
    dist += D[t];
    return D[t] < inf;
}

long long Fuck::costflow() {
    flow = dist = 0;
    cost = 0;
    while (modlabel()) {
        int tmpf;
        do {
            memset(vis , 0, sizeof(vis));
            tmpf = aug(s, INT_MAX);
            flow += tmpf;
        } while (tmpf);
    }
    for (Edge *i = e[s]; i; i = i->next) {
        cost += node_cost[i->t] + flow_cost[i->U - i->u].cost;
    }
    if (flow != need_flow)
        cost = infll;
    return cost;
}


Particle::Particle(int length): v(length, 0), v_best(length, 0), vp(length, 0), cost_best(infll), cost(infll) {}

Particle::Particle(int length, vector<int> & vi, Fuck * & fuck): v(length, 0), v_best(length, 0), vp(length, 0) {
    int size = vi.size();
    for (int i = 0; i < size; ++i) {
        v_best[vi[i]] = v[vi[i]] = 1;
    }
    fuck->add_server(vi);
    cost = cost_best = fuck->costflow();
}
bool cmp(const Particle & p1, const Particle & p2) {
    return p1.cost == p2.cost ? p1.cost_best < p2.cost_best : p1.cost < p2.cost;
}

XJBS::XJBS(Fuck & fk) {
    fuck = &fk;
    l = fuck->node_num;
}

inline void XJBS::decode(vector<double> & vd, vector<int> & vi) {
    vi.clear();
    for (int i = 0; i < l; ++i) {
        if (vd[i] > 0.5)
            vi.push_back(i);
    }
}

inline void XJBS::addone(vector<int> & v) {
    p.emplace_back(l, v, fuck);
}

inline void XJBS::get_best(vector<int> & server, long long & cost) {
    decode(gbest.v_best, server);
    cost = gbest.cost_best;
}

inline void XJBS::GA_cross(Particle & s1, Particle & s2) {
    //clock_t t1 = clock();
    int r1 = rand() % l, r2 = rand() % l;
    if (r2 < r1)
        swap(r1, r2);
    while (r1 < r2) {
        swap(s1.v[r1], s2.v[r1]);
        ++r1;
    }
    //cout << "GA_cross:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}


inline void XJBS::OBMA(Particle & s) {
    //clock_t t1 = clock();
    int r1, r2;
    do {
        r1 = rand() % l;
    } while (s.v[r1] > 0.5);
    do {
        r2 = rand() % l;
    } while (s.v[r2] < 0.5);
    swap(s.v[r1], s.v[r2]);
    
    vector<int> server, unserver;
    for (int i = 0; i < l; ++i) {
        if (s.v[i] > 0.5) {
            server.push_back(i);
        } else {
            unserver.push_back(i);
        }
    }
    int server_size = server.size(), unserver_size = unserver.size(), server_no, unserver_no;
    knuth_shuffle(server);
    knuth_shuffle(unserver);
    if (server_size && unserver_size) {
        swap(s.v[server[0]], s.v[unserver[0]]);
    }
    for (int i = 0, server_index = 0, unserver_index = 0; i < 1 && server_index < server_size && unserver_index < unserver_size; ++i, ++server_index, ++unserver_index) {
        server_no = server[server_index];
        unserver_no = unserver[unserver_index];
        swap(s.v[server_no], s.v[unserver_no]);
        swap(server[server_index], unserver[unserver_index]);
    }
    //cout << "OBMA:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}

inline void XJBS::PSO_update(Particle & s) {
    //clock_t t1 = clock();
    for (int i = 0; i < l; ++i) {
        s.vp[i] = PSO_w * s.vp[i] + PSO_c1 * rand() / RAND_MAX * (s.v_best[i] - s.v[i]) + PSO_c2 * rand() / RAND_MAX * (gbest.v_best[i] - s.v[i]); 
        s.v[i] = (1 / (1 + exp(100*(0.5-(s.v[i] + s.vp[i])))));
    }
    //cout << "PSO:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}

inline void XJBS::updateone(Particle & s) {
    vector<int> v;
    decode(s.v, v);
    fuck->add_server(v);
    s.cost = fuck->costflow();
    if (s.cost < s.cost_best) {
        s.v_best = s.v;
        s.cost_best = s.cost;
        if (s.cost_best < gbest.cost_best) {
            gbest.v_best = s.v_best;
            gbest.cost_best = s.cost_best;
            cnt = 0;
        }
    }
}

inline void XJBS::reproduction() {
    for (int i = 0; i < max_p_size; ++i)
        p.push_back(p[i]);
    max_p_size >>= 1;
}

void XJBS::run1() {
    for (int i = 0; i < max_p_size; ++i) {
        OBMA(p[i]);
        updateone(p[i]);
        PSO_update(p[i]);
    }
    if (++cnt > 100) {
        run_second >>= 1;
        cnt = 0;
    }
    //cout << "run1 " << gbest.cost_best << endl;
}

void XJBS::run2() {
    int i = 0;
    int j = max_p_size - 1;
    sort(p.begin(), p.end(), cmp);
    for (; i < j; ++i, --j) {
        OBMA(p[i]);
        PSO_update(p[j]);
        GA_cross(p[i], p[j]);
        updateone(p[i]);
        updateone(p[j]);
    }
    if (++cnt > 100) {
        run_second >>= 1;
        cnt = 0;
    }
    //cout << "run2 " << gbest.cost_best << endl;
}

void XJBS::initial() {
    PSO_c1 = c1;
    PSO_c2 = c2;
    PSO_w = w;
    cnt = 0;
    sort(p.begin(), p.end(), cmp);
    gbest = p[0];
    int best_size = 0;;
    for (int i = 0; i < l; ++i)
        best_size += gbest.v_best[i] > 0.5 ? 1 : 0;
    if (best_size > 150)
        max_p_size = 8;
    else if (best_size > 100)
        max_p_size = 16;
    else 
        max_p_size = 20;
    best_size *= 0.7;
    int limit_size = min(max_p_size >> 1, (int)p.size());
    p.resize(limit_size);
    vector<int> v;
    for (int i = limit_size; i < max_p_size; ++i) {
        fuck->kmeans(best_size, v);
        addone(v);
    }
}

