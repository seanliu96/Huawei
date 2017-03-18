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

using namespace std;

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    char * topo_file;
    LiuXin liuxin;
    liuxin.readtopo(topo, line_num);
    liuxin.spfa();
    JinTao jintao(liuxin);
    HGAPSO hgapso(jintao, pm, pc, c1, c2, w);
    vector<int> best_server;
    long long best_cost = infll;
    vector<vector<int> > node;
    vector<int> flow;
    int best_index = liuxin.customer_num;
    int block_size = (int)(sqrt(liuxin.customer_num) + 0.1);
    
    for (int i = 1; i <= liuxin.customer_num; i += block_size) {
        vector<int> server1 = liuxin.kmeans(i), server2;
        jintao.recover();
        jintao.add_server(server1);
        long long best_now = jintao.costflow() + i * (long long)liuxin.server_cost;
        if (best_now < best_cost) {
            best_cost = best_now;
            best_server = server1;
            best_index = i;
        }
        for (int j = 1; j < kmean_times ; ++j) {
            server2 = liuxin.kmeans(i);
            bool same = true;
            for (int k = 0; k < i; ++k) {
                if (server1[k] != server2[k]) {
                    same = false;
                    break;
                }
            }
            if (same)
                continue;
            jintao.recover();
            jintao.add_server(server2);
            long long cost = jintao.costflow() + i * (long long)liuxin.server_cost;
            if (cost < best_now) {
                server1.swap(server2);
                best_now = cost;
                if (best_now < best_cost) {
                    best_cost = best_now;
                    best_server = server1;
                    best_index = i;
                }
            }
        }
    }
    block_size >>= 1;
    for (int i = max(best_index - block_size, 1); i <= min(best_index + block_size, liuxin.customer_num); i++) {
        vector<int> server1 = liuxin.kmeans(i), server2;
        hgapso.addone(server1);
        //jintao.recover();
        //jintao.add_server(server1);
        /*
        long long best_now = jintao.costflow() + i * (long long)liuxin.server_cost;
        if (best_now < best_cost) {
            best_cost = best_now;
            best_server = server1;
            best_index = i;
        }
        */
        for (int j = 1; j < kmean_times ; ++j) {
            server2 = liuxin.kmeans(i);
            bool same = true;
            for (int k = 0; k < i; ++k) {
                if (server1[k] != server2[k]) {
                    same = false;
                    break;
                }
            }
            if (same)
                continue;
            hgapso.addone(server2);
            server1.swap(server2);
            /*
            jintao.recover();
            jintao.add_server(server2);
            */
            /*
            long long cost = jintao.costflow() + i * (long long)liuxin.server_cost;
            if (cost < best_now) {
                server1 = server2;
                hgapso.addone(server1);
                best_now = cost;
                if (best_now < best_cost) {
                    best_cost = best_now;
                    best_server = server1;
                    best_index = i;
                }
            }
            */
        }
    }
    int hgapso_times = hgapso.first_run() << 3;
    while (hgapso.run() < hgapso_times && (double)clock() / CLOCKS_PER_SEC < last_second);
    best_server = hgapso.get_best();
    jintao.recover();
    jintao.add_server(best_server);
    jintao.costflow();
    jintao.print_flow(node, flow);
    topo_file = new char[node.size() * MAX_V * 4];
    topo_file[0] = '\0';
    char tmp[MAX_V * 4];
    sprintf(tmp, "%d\n\n", (int)node.size());
    strcat(topo_file, tmp);
    for (unsigned int i = 0; i < node.size(); i++) {
        for (unsigned int j = 0; j< node[i].size(); j++) {
            if (j == node[i].size()-1)
                node[i][j] -= liuxin.node_num;
            sprintf(tmp, "%d ", node[i][j]);
            strcat(topo_file, tmp);
        }
        sprintf(tmp, "%d\n", (int)flow[i]);
        strcat(topo_file, tmp);
    }
    write_result(topo_file, filename);
    delete []topo_file;
}

Particle::Particle(int length): v(length, 0.0), v_best(length, 0.0), vp(length, 0.0), cost_best(infll), cost(infll) {}
/*
bool operator== (const Particle & p1, const Particle & p2) {
    return p1.cost_best == p2.cost_best && p1.cost == p2.cost;
}

bool operator< (const Particle & p1, const Particle & p2) {
    return p1.cost_best == p2.cost_best ? p1.cost < p2.cost : p1.cost_best < p2.cost_best;
}

bool operator> (const Particle & p1, const Particle & p2) {
    return p1.cost_best == p2.cost_best ? p1.cost > p2.cost : p1.cost_best > p2.cost_best;
}

bool operator<= (const Particle & p1, const Particle & p2) {
    return p1 < p2 || p1 == p2;
}

bool operator>= (const Particle & p1, const Particle & p2) {
    return p1 > p2 || p1 == p2;
}
*/


bool operator== (const Particle & p1, const Particle & p2) {
    return p1.cost_best == p2.cost_best && p1.cost == p2.cost;
}

bool operator< (const Particle & p1, const Particle & p2) {
    return p1.cost == p2.cost ? p1.cost_best < p2.cost_best : p1.cost < p2.cost;
}

bool operator> (const Particle & p1, const Particle & p2) {
    return p1.cost == p2.cost ? p1.cost_best > p2.cost_best : p1.cost > p2.cost;
}

bool operator<= (const Particle & p1, const Particle & p2) {
    return p1 < p2 || p1 == p2;
}

bool operator>= (const Particle & p1, const Particle & p2) {
    return p1 > p2 || p1 == p2;
}


template <class T>
void knuth_shuffle(vector<T> & v) {
    int i = v.size() - 1, j = 0;
    while (i >= 0) {
        j = rand() % (i + 1);
        swap(v[i], v[j]);
        --i;
    }
}

void LiuXin::readtopo(char * topo[MAX_EDGE_NUM], int line_num) { 
    int line = 0;
    int u, v, c, w;
    if (line < line_num)
        sscanf(topo[line], "%d %d %d", &node_num, &edge_num, &customer_num);
    graph.resize(node_num, vector<EdgeInfo>());
    d.resize(node_num, vector<int>(node_num, inf));
    customer_nodes.resize(customer_num);
    line+=2;
    sscanf(topo[line], "%d", &server_cost);
    line+=2;
    for (int i = 0; i < edge_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d %d", &u, &v, &w, &c);
        graph[u].push_back((EdgeInfo){v, w, c});
        graph[v].push_back((EdgeInfo){u, w, c});
    }
    ++line;
    need_flow = 0;
    for (int i = 0; i < customer_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d", &u, &v, &w);
        customer_nodes[i] = (CustomerNodeInfo){u, v, w};
        need_flow += w;
    }
}

void LiuXin::spfa() {
    for(int s = 0; s < node_num; ++s) {
        d[s][s] = 0;
        queue<int> Q;
        memset(vis,0,sizeof(vis));
        vis[s] = true;
        Q.push(s);
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();
            for (unsigned int i = 0; i < graph[u].size(); ++i) {
                int v = graph[u][i].v;
                int dis =  d[s][u] + graph[u][i].c;
                if (dis < d[s][v])
                { 
                    d[s][v] = dis;
                    if (!vis[v])
                    {
                        vis[v] = true;
                        Q.push(v);
                    }
                }
            }
            vis[u] = false;
        }
    }
}

vector<int> LiuXin::kmeans(int k) {
    srand(time(0));
    vector<int> clusters(k);
    int label[customer_num], min_dist, min_index;
    knuth_shuffle(customer_nodes);
    for (int i = 0; i < k; ++i) {
        clusters[i] = customer_nodes[i].v;
    }
    bool update = true;
    while (update) {
        update = false;
        for (int i = 0; i < customer_num; ++i) {
            min_dist= inf;
            min_index = 0;
            for (int j = 0; j < k; ++j) {
                if(d[clusters[j]][customer_nodes[i].v] < min_dist) {
                    min_dist = d[clusters[j]][customer_nodes[i].v];
                    min_index = j;
                }
            }
            if (label[i] != min_index) {
                update = true;
                label[i] = min_index;
            }
        }
        for (int j = 0; j < k; ++j) {
            min_dist = inf;
            min_index = -1;
            for (int l = 0; l < node_num; ++l) {
                int dist = 0;
                for (int i = 0; i < customer_num; ++i) {
                    dist += label[i] == j ? d[l][customer_nodes[i].v] : 0;
                }
                if (dist < min_dist) {
                    min_index = l;
                    min_dist = dist;
                }
            }
            clusters[j] = min_index;
        }
    }
    sort(clusters.begin(), clusters.end());
    return clusters;
}

JinTao::JinTao(LiuXin & lx) {
    liuxin = &lx;
    recover();
}

void JinTao::recover() {
    psz = 0;
    memset(e, 0, sizeof(e));
    for (unsigned int u = 0; u < liuxin->graph.size(); ++u) {
        for (unsigned int i = 0; i < liuxin->graph[u].size(); ++i) {
            add_edge(u, liuxin->graph[u][i].v, liuxin->graph[u][i].w, liuxin->graph[u][i].c);
        }
    }
    s = liuxin->node_num + liuxin->customer_num;
    t = s + 1;
    for (int i = 0; i < liuxin->customer_num; ++i) {
        add_edge(liuxin->customer_nodes[i].v, liuxin->customer_nodes[i].u + liuxin->node_num, liuxin->customer_nodes[i].w, 0);
        add_edge(liuxin->customer_nodes[i].u + liuxin->node_num, t, liuxin->customer_nodes[i].w, 0);
    }
}

void JinTao::add_edge(int u, int v, int w, int c) {
    Edge *e1 = epool + psz++, *e2 = epool + psz ++;
    *e1 = (Edge){v, w, c, w, e[u], e2}, e[u] = e1;
    *e2 = (Edge){u, 0, -c, 0, e[v], e1}, e[v] = e2;
}


void JinTao::add_server(vector<int> &Q) {
    for (unsigned int i = 0; i < Q.size(); i++) {
        add_edge(s, Q[i], inf, 0);
    }
}

void JinTao::print_flow(vector<vector<int> > & node, vector<int> &flow) {
    while (true) {
        vector<int> Tmp;
        int u = s;
        int S = inf;
        while (u != t) {
            bool flag=false;
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    S = min(S, i->U - i->u);
                    u = v;
                    flag=true;
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

int JinTao::aug(int u, int m) {
    if (u == t)
        return cost += (long long)dist * m, m;
    int d = m;
    vis[u] = true;
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

bool JinTao::modlabel() {
    deque <int > q;
    memset(vis , 0, sizeof(vis));
    memset(d, 0x3f , sizeof(d));
    q.push_back(s);
    d[s] = 0;
    vis[s] = true;
    while (!q.empty ()) {
        int u = q.front ();
        q.pop_front ();
        vis[u] = false;
        for (Edge *i = e[u]; i; i = i->next) {
            int v = i->t;
            if (i->u && d[u] + i->c < d[v]) {
                d[v] = d[u] + i->c;
                if (vis[v])
                    continue;
                vis[v] = true;
                if (q.size () && d[v] < d[q[0]])
                    q.push_front(v);
                else
                    q.push_back(v);
            }
        }
    }
    for (Edge *i = epool; i < epool + psz; ++i) {
        i->c -= d[i->t] - d[i->pair->t];
    }
    dist += d[t];
    return d[t] < inf;
}

long long JinTao::costflow() {
    flow = dist = 0;
    cost = 0;
    while (modlabel ()) {
        int tmpf;
        do {
            memset(vis , 0, sizeof(vis));
            tmpf = aug(s, INT_MAX);
            flow += tmpf;
        } while (tmpf);
    }
    if (flow != liuxin->need_flow)
        cost = infll;
    return cost;
}

HGAPSO::HGAPSO(JinTao & jt, double pm, double pc, double c1, double c2, double w) {
    GA_pm = pm;
    GA_pc = pc;
    PSO_c1 = c1;
    PSO_c2 = c2;
    PSO_w = w;
    jintao = &jt;
    l = jintao->liuxin->node_num;
    gbest = Particle(l);
    unchanged_times = 0;
}

Particle HGAPSO::encode(vector<int> & v) {
    Particle res = Particle(l);
    for (unsigned int i = 0; i < v.size(); ++i) {
        res.v[v[i]] = 1.0;
    }
    return res;
}

vector<int> HGAPSO::decode(vector<double> & v) {
    vector<int> res;
    for (int i = 0; i < l; ++i) {
        if (v[i] > 0.5)
            res.push_back(i);
    }
    return res;
}

void HGAPSO::initial() {
    p.clear();
    unchanged_times = 0;
    gbest = Particle(l);
}

void HGAPSO::addone(vector<int> & v) {
    p.push_back(encode(v));
}

vector<int> HGAPSO::get_best() {
    return decode(gbest.v_best);
}

void HGAPSO::GA_cross(Particle & s1, Particle & s2) {
    /*
    for (int i = 0; i < l; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < pc) {
            swap(s1.v[i], s2.v[i]);
        }   
    }
    */
    
    double r = (double)rand() / RAND_MAX;
    if (r > GA_pc)
        return;
    int r1 = rand() % l, r2 = rand() % l;
    if (r1 > r2) 
        swap(r1, r2);
    while (r1 < r2) {
        swap(s1.v[r1], s2.v[r1]);
        ++r1;
    }
    
    /*
    int r1 = rand() % l, r2;
    while (r1--) {
        r2 = rand() % l;
        swap(s1.v[r2], s2.v[r2]);
    }
    */
}

void HGAPSO::GA_mutation(Particle & s) {
    
    for (int i = 0; i < l; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < GA_pm) {
            s.v[i] += ((double)rand() / RAND_MAX) * 2 - 1;
            s.v[i] = max(0.0, s.v[i]);
            s.v[i] = min(s.v[i], 1.0);
        }   
    }
    /*
    if ((double)rand() / RAND_MAX > GA_pm)
        return;
    int r1 = rand() % (gbest.v_best.size()), r2;
    while (r1--) {
        r2 = rand() % l;
        s.v[r2] += (double)rand() / RAND_MAX * 2 - 1;
        s.v[r2] = max(0.0, s.v[r2]);
        s.v[r2] = min(s.v[r2], 1.0);
    }
    */
}

void HGAPSO::PSO_update(Particle & s) {
    for (int i = 0; i < l; ++i) {
        s.vp[i] = PSO_w * s.vp[i] + PSO_c1 * rand() / RAND_MAX * (s.v_best[i] - s.v[i]) + PSO_c2 * rand() / RAND_MAX * (gbest.v_best[i] - s.v[i]); 
        s.v[i] = s.v[i] + s.vp[i];
        s.v[i] = max(0.0, s.v[i]);
        s.v[i] = min(s.v[i], 1.0);
    }
}

int HGAPSO::run() {
    int r1, r2, i = 0, k = p.size(), j = k >> 1;
    for (int i = 0; i < k; ++i) {
        jintao->recover();
        vector<int> v = decode(p[i].v);
        jintao->add_server(v);
        p[i].cost = (long long)jintao->costflow() + v.size() * (long long)jintao->liuxin->server_cost;
        if (p[i].cost < p[i].cost_best) {
            p[i].v_best = p[i].v;
            p[i].cost_best = p[i].cost;
            if (p[i].cost < gbest.cost_best) {
                gbest.v_best = p[i].v_best;
                gbest.cost_best = p[i].cost_best;
                unchanged_times = 0;
            }
        }
    }
    for (i = 0; i < k; ++i) {
        PSO_update(p[i]);
    }
    sort(p.begin(), p.end());
    for (i = j; i < k; ++i) {
        r1 = rand() % j;
        r2 = rand() % j;
        p[i] = p[r1] < p[r2] ? p[r1] : p[r2];
        GA_cross(p[i - 1], p[i]);
        GA_mutation(p[i]);
    }
    /*
    for (i = 1; i <= 10; i++) {
        r1 = rand() % (k - j);
        r2 = rand() % (k - j);
        GA_cross(p[r1 + j], p[r2 + j]);
    }
    */
    ++unchanged_times;
    cout << unchanged_times << endl;
    return unchanged_times;
}

int HGAPSO::first_run() {
    run();
    unchanged_times = 0;
    int best_size = get_best().size();
    //int new_size = best_size << 2;
    int new_size = MAX_P_SIZE;
    vector<int> v;
    for (unsigned int  i = p.size(); i < new_size; ++i) {
        v = jintao->liuxin->kmeans(best_size);
        addone(v);
    }
    p.resize(new_size);
    return new_size;
}
