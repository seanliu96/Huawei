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

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    char * topo_file;
    srand(time(0));
    Fuck fuck;
    fuck.readtopo(topo, line_num);
    fuck.spfa();
    HGAPSO hgapso(fuck);
    vector<int> server, best_server;
    vector<vector<int> > node;
    vector<int> flow;
    int best_index = fuck.customer_num;
    int kmean_times = 2;
    int block_size = (int)(sqrt(fuck.customer_num) + 0.1);
    double last_second = 90 - 1;
    fuck.kmeans(1, server);
    hgapso.addone(server);
    fuck.kmeans(best_index, best_server);
    long long best_cost = best_index * fuck.server_cost;
    for (int i = best_index - 1; i > 1; i -= block_size) {
        for (int j = 0; j < kmean_times ; ++j) {
            fuck.kmeans(i, server);
            fuck.add_server(server);
            long long cost = fuck.costflow() + i * fuck.server_cost;
            if (cost < best_cost) {
                best_cost = cost;
                best_index = i;
                best_server.swap(server);
            }
        }
    }
    hgapso.addone(best_server);
    block_size = (block_size >> 1);
    int min_index = max(best_index - block_size, 1);
    int max_index = min(best_index + block_size, fuck.customer_num);
    int max_p_size = 128 / log(best_index << 3);
    max_p_size += (max_p_size & 1);
    for (int i = min_index; i <= max_index; ++i) {
        for (int j = 0; j < kmean_times; ++j) {
            fuck.kmeans(i, server);
            hgapso.addone(server);
        }
    }
    last_second -= hgapso.initial(max_p_size);
    while ((double)clock() / CLOCKS_PER_SEC < last_second)
        hgapso.run();
    hgapso.get_best(best_server);
    fuck.add_server(best_server);
    fuck.costflow();
    fuck.print_flow(node, flow);
    int node_size = node.size();
    topo_file = new char[node_size * MAX_V * 5];
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
        sprintf(tmp, "%d\n", (int)flow[i]);
        strcat(line, tmp);
        strcat(topo_file, line);
    }
    write_result(topo_file, filename);
    delete []topo_file;
}

Particle::Particle(int length): v(length, 0.0), v_best(length, 0.0), vp(length, 0.0), cost_best(infll), cost(infll) {}

Particle::Particle(int length, vector<int> & vi, Fuck * & fuck): v(length, 0.0), v_best(length, 0.0), vp(length, 0.0), cost_best(infll), cost(infll) {
    int size = vi.size();
    for (int i = 0; i < size; ++i) {
        v_best[vi[i]] = v[vi[i]] = 1.0;
    }
    fuck->add_server(vi);
    cost = cost_best = fuck->costflow() + size * fuck->server_cost;
}
bool cmp(const Particle & p1, const Particle & p2) {
    return p1.cost == p2.cost ? p1.cost_best < p2.cost_best : p1.cost < p2.cost;
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
    sscanf(topo[line], "%lld", &server_cost);
    line += 2;
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
        vis[s] = true;
        q.push_back(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            vis[u] = false;
            for (unsigned int i = 0; i < graph[u].size(); ++i) {
                int v = graph[u][i].v;
                int dis =  d[s][u] + graph[u][i].c;
                if (dis < d[s][v])
                { 
                    d[s][v] = dis;
                    if (!vis[v]) {
                        vis[v] = true;
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
    int label[MAX_V];
    memset(label, -1, sizeof(label));
    vector<vector<int> > kmean_node(k);
    int min_dist, min_index;
    knuth_shuffle(customer_nodes);
    for (int i = 0; i < k; ++i) {
        clusters[i] = customer_nodes[i].v;
    }
    while (true) {
        for (int i = 0; i < k; ++i)
            kmean_node[i].clear();
        bool update = false;
        for (int i = 0; i < customer_num; ++i) {
            min_dist = inf;
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
            kmean_node[label[i]].push_back(i);
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
        add_edge(s, v[i], inf, 0);
    }
}

void Fuck::add_edge(int u, int v, int w, int c) {
    Edge *e1 = epool + psz++, *e2 = epool + psz++;
    *e1 = (Edge){v, w, c, w, c, e[u], e2}, e[u] = e1;
    *e2 = (Edge){u, 0, -c, 0, -c, e[v], e1}, e[v] = e2;
}


void Fuck::print_flow(vector<vector<int> > & node, vector<int> &flow) {
    node.clear();
    flow.clear();
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
                    flag = true;
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

bool Fuck::modlabel() {
    deque <int> q;
    memset(vis , 0, sizeof(vis));
    memset(D, 0x3f , sizeof(D));
    q.push_back(s);
    D[s] = 0;
    vis[s] = true;
    while (!q.empty ()) {
        int u = q.front ();
        q.pop_front ();
        for (Edge *i = e[u]; i; i = i->next) {
            int v = i->t;
            int dis = D[u] + i->c;
            if (i->u && dis < D[v]) {
                D[v] = dis;
                if (!vis[v]) {
                    vis[v] = true;
                    if (q.size () && D[v] < D[q[0]])
                        q.push_front(v);
                    else
                        q.push_back(v);
                }
            }
        }
        vis[u] = false;
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
    if (flow != need_flow)
        cost = infll;
    return cost;
}

HGAPSO::HGAPSO(Fuck & fk) {
    PSO_c1 = c1;
    PSO_c2 = c2;
    PSO_w = w;
    fuck = &fk;
    l = fuck->node_num;
    gbest = Particle(l);
    H.resize(l, 0);
    iter = 0;
}

void HGAPSO::decode(vector<double> & vd, vector<int> & vi) {
    vi.clear();
    for (int i = 0; i < l; ++i) {
        if (vd[i] > 0.5)
            vi.push_back(i);
    }
}

void HGAPSO::addone(vector<int> & v) {
    p.emplace_back(l, v, fuck);
}

void HGAPSO::get_best(vector<int> & server) {
    decode(gbest.v_best, server);
}

void HGAPSO::cross(Particle & s1, Particle & s2) {
    //clock_t t1 = clock();
    vector<double> v1(l, 0), v2(l, 0);
    int m = 0, l_2 = l >> 1;
    vector<int> node1, node2;
    long long cost;
    for (int i = 0; i < l; ++i) {
        v1[i] = (s1.v[i] + s2.v[i]) * 0.5;
        if (v1[i] > 0.5)
            node1.push_back(i);
        else
            node2.push_back(i);
    }
    m = node1.size();
    knuth_shuffle(node1);
    knuth_shuffle(node2);
    if (m < l_2) {
        for (int i = 0; i < m; ++i)
            v2[node2[i]] = 1 - v1[node2[i]];
    } else if (m == l_2) {
        for (int i = 0; i < l; ++i)
            v2[i] = 1 - v1[i];
    } else {
        int j = node2.size();
        for (int i = 0; i < j; ++i)
            v2[node2[i]] = 1 - v1[node2[i]];
        j = (m << 1) - j;
        for (int i = 0; i < j; ++i)
            v2[node1[i]] = v1[node1[i]];
    }
    decode(v1, node1);
    fuck->add_server(node1);
    cost = fuck->costflow() + node1.size() * fuck->server_cost;
    if (cost < s1.cost) {
        s1.cost = cost;
        s1.v = v1;
        if (cost < s1.cost_best) {
            s1.cost_best = cost;
            s1.v_best = v1;
            if (cost < gbest.cost_best) {
                gbest.cost_best = cost;
                gbest.v_best = v1;
                PSO_c1 = c1;
                PSO_c2 = c2;
                PSO_w = w;
            }
        }
    }
    decode(v2, node2);
    fuck->add_server(node2);
    cost = fuck->costflow() + node1.size() * fuck->server_cost;
    if (cost < s2.cost) {
        s2.cost = cost;
        s2.v = v2;
        if (cost < s2.cost_best) {
            s2.cost_best = cost;
            s2.v_best = v2;
            if (cost < gbest.cost_best) {
                gbest.cost_best = cost;
                gbest.v_best = v2;
                PSO_c1 = c1;
                PSO_c2 = c2;
                PSO_w = w;
            }
        }
    }
    //cout << "cross:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}

void HGAPSO::OBMA(Particle s) {
    //clock_t t1 = clock();
    vector<double> v(s.v);
    vector<int> server, unserver;
    for (int i = 0; i < l; ++i) {
        if (v[i] > 0.5)
            server.push_back(i);
        else
            unserver.push_back(i);
    }
    long long best_cost = s.cost, cost;
    int server_index, unserver_index, server_size = server.size(), unserver_size = unserver.size(), server_no, unserver_no;
    for (int k = 0; k < 5; ++k) {
        ++iter;
        int T_iter = 15 * ((iter >> 8) + 1);
        do {
            server_index = rand() % server_size;
            server_no = server[server_index];
        } while (iter < H[server_no]);
        H[server_no] = T_iter;
        do {
            unserver_index = rand() % unserver_size;
            unserver_no = unserver[unserver_index];
        } while (iter < H[unserver_no]);
        H[unserver_no] = 0.7 * T_iter;
        swap(v[server_no], v[unserver_no]);
        swap(server[server_index], unserver[unserver_index]);
        fuck->add_server(server);
        cost = fuck->costflow() + server.size() * fuck->server_cost;
        if (cost < best_cost) {
            s.v = v;
            s.cost = cost;
            if (cost < s.cost_best) {
                s.v_best = v;
                s.cost_best = cost;
            }
            p.push_back(s);
        }
    }
    if (s.cost_best < gbest.cost_best) {
        gbest.v_best = s.v_best;
        gbest.cost_best = s.cost_best;
        PSO_c1 = c1;
        PSO_c2 = c2;
        PSO_w = w;
    }
    //cout << "OBMA:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}

void HGAPSO::PSO_update(Particle & s) {
    //clock_t t1 = clock();
    vector<int> v;
    long long cost;
    for (int i = 0; i < l; ++i) {
        s.vp[i] = PSO_w * s.vp[i] + PSO_c1 * rand() / RAND_MAX * (s.v_best[i] - s.v[i]) + PSO_c2 * rand() / RAND_MAX * (gbest.v_best[i] - s.v[i]); 
        //s.v[i] = 1 / (1 + exp(0.5 - s.v[i] - s.vp[i]));
        s.v[i] = s.v[i] + s.vp[i];
        s.v[i] = max(0.0, s.v[i]);
        s.v[i] = min(s.v[i], 1.0);
    }
    decode(s.v, v);
    fuck->add_server(v);
    cost = fuck->costflow() + v.size() * fuck->server_cost;
    s.cost = cost;
    if (cost < s.cost_best) {
        s.v_best = s.v;
        s.cost_best = s.cost;
        if (cost < gbest.cost_best) {
            gbest.v_best = s.v;
            gbest.cost_best = cost;
            PSO_c1 = c1;
            PSO_c2 = c2;
            PSO_w = w;
        }
    }
    //cout << "PSO:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}

void HGAPSO::run() {
    int r1, r2, i, k = p.size(), j = k * 3 >> 2;
    vector<int > v;
    
    for (i = 0; i < k; ++i) {
        PSO_update(p[i]);
        OBMA(p[i]);
    }
    sort(p.begin(), p.end(), cmp);
    p.resize(max_p_size);
    /* 
    sort(p.begin(), p.end(), cmp);
    for (i = k - 1; i >=j; i -= 2) {
        r1 = rand() % j;
        do {
            r2 = rand() % j;
        } while (r1 == r2);
        p[i-1] = cmp(p[r1], p[r2]) ? p[r1] : p[r2];
        r1 = rand() % j;
        do {
            r2 = rand() % j;
        } while (r1 == r2);
        p[i] = cmp(p[r1], p[r2]) ? p[r1] : p[r2];
        cross(p[i-1], p[i]);
    }
    for (i = 0; i < k; ++i)
        OBMA(p[i]);
    */
    PSO_c1 *= alpha;
    PSO_c2 *= alpha;
    PSO_w *= alpha;
    //cout << gbest.cost_best << endl;
}

double HGAPSO::initial(int size) {
    max_p_size = size;
    int p_size = p.size(), limit_size = max_p_size * 3 >> 2;
    vector<int> v;
    sort(p.begin(), p.end(), cmp);
    gbest = p[0];
    decode(gbest.v_best, v);
    int best_size = v.size();
    if (p_size < limit_size) {
        for (int i = p_size; i < max_p_size; ++i) {
            fuck->kmeans(best_size, v);
            addone(v);
        }
    } else {
        p.resize(limit_size);
        for (int i = limit_size; i < max_p_size; ++i) {
            fuck->kmeans(best_size, v);
            addone(v);
        }
    }
    clock_t t1 = clock();
    run();
    clock_t t2 = clock();
    return double(t2 - t1) / CLOCKS_PER_SEC;
}
