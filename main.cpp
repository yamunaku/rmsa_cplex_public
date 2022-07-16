#include <ilcplex/ilocplex.h>
#include <vector>
#include <random>
#include <queue>
#include <cmath>
#include <algorithm>

using namespace std;

struct Edge {
    int from = -1;
    int to = -1;
    int cost = 0;

    Edge() {}
    Edge(int _from, int _to, int _cost) : from(_from), to(_to), cost(_cost) {}

    int next(int v) const {
        if (v != from && v != to) return -1;
        return from + to - v;
    }
};

struct Graph {
private:
    int vertices_num = 0;
    vector<Edge> edges;

public:
    Graph() {}
    Graph(int _vertices_num) : vertices_num(_vertices_num), edges(0) {}
    Graph(int _vertices_num, vector<Edge> _edges) : vertices_num(_vertices_num), edges(_edges) {}

    void set_vertices_num(int _vertices_num) {
        vertices_num = _vertices_num;
    }

    int get_vertices_num() const {
        return vertices_num;
    }

    void set_edges(vector<Edge> _edges) {
        edges = _edges;
    }

    void add_edge(Edge edge) {
        edges.push_back(edge);
    }

    vector<Edge> get_edges() const {
        return edges;
    }

    pair<int, vector<int>> dijkstra(const vector<vector<int>> &edges_lists, const int s, const int t, const vector<int> &vertices_usage, const vector<int> &edges_usage) const {
        const int INF = 1000000000;
        int n = vertices_num;
        int m = edges.size();
        vector<int> can_use_vertex(n, true);
        vector<int> can_use_edge(m, true);
        for (int x : vertices_usage) can_use_vertex[x] = false;
        for (int x : edges_usage) can_use_edge[x] = false;
        vector<int> dist(n, INF), par(n, -1);
        priority_queue<pair<int, pair<int, int>>> heap;
        heap.push({ 0, {s, -1} });
        while (!heap.empty()) {
            auto p = heap.top();
            heap.pop();
            int now = p.second.first;
            if (dist[now] != INF) continue;
            dist[now] = -p.first;
            par[now] = p.second.second;
            for (auto e_id : edges_lists[now]) {
                if (!can_use_edge[e_id]) continue;
                auto e = edges[e_id];
                int nx = e.next(now);
                if (!can_use_vertex[nx]) continue;
                if (dist[nx] != INF) continue;
                heap.push({ -dist[now] - e.cost, {nx, e_id} });
            }
        }

        if (dist[t] == INF) {
            return { -1, {} };
        }

        pair<int, vector<int>> ret;
        ret.first = dist[t];
        int now = t;
        while (par[now] != -1) {
            int e_id = par[now];
            ret.second.push_back(e_id);
            now = edges[e_id].next(now);
        }
        reverse(ret.second.begin(), ret.second.end());

        return ret;
    }

    vector<vector<int>> yen_algorithm(const vector<vector<int>>& edges_lists, const int s, const int t, const int k_max) const {
        int n = vertices_num;
        vector<vector<int>> route(k_max);

        priority_queue<pair<int, pair<pair<int,int>, vector<int>>>> heap;
        vector<int> par(k_max), branch(k_max);
        auto path_0 = dijkstra(edges_lists, s, t, {}, {});
        heap.push({ -path_0.first, {{-1, 0}, path_0.second} });

        for (int k = 0; k < k_max; k++) {
            if (heap.empty()) break;
            auto p_k = heap.top();
            heap.pop();
            route[k] = p_k.second.second;
            par[k] = p_k.second.first.first;
            branch[k] = p_k.second.first.second;

            if (k == k_max - 1) break;
            vector<int> vertices_usage = {};
            int now = s;
            int cost_sum = 0;
            vector<int> path;
            for (int i = 0; i < (int)route[k].size(); i++) {
                int e_id = route[k][i];
                vertices_usage.push_back(now);

                if (i >= branch[k]) {
                    vector<int> edges_usage = { route[k][i] };
                    if (i == branch[k]) {
                        int now_k = par[k];
                        while (now_k != -1) {
                            edges_usage.push_back(route[now_k][i]);
                            if (branch[now_k] != branch[k]) break;
                            now_k = par[now_k];
                        }
                    }
                    auto ret = dijkstra(edges_lists, now, t, vertices_usage, edges_usage);
                    if (ret.first != -1) {
                        vector<int> p = path;
                        p.insert(p.end(), ret.second.begin(), ret.second.end());
                        heap.push({ -cost_sum - ret.first, {{k, i}, p} });
                    }
                }

                auto e = edges[e_id];
                path.push_back(e_id);
                cost_sum += e.cost;
                now = e.next(now);
            }
        }

        return route;
    }

    vector<vector<int>> make_edges_lists() const {
        vector<vector<int>> edges_lists(vertices_num);
        for (int i = 0; i < edges.size(); i++) {
            edges_lists[edges[i].from].push_back(i);
            edges_lists[edges[i].to].push_back(i);
        }
        return edges_lists;
    }

    vector<vector<vector<vector<int>>>> calc_route(const int k_max) const {
        int n = vertices_num;
        const auto edges_lists = make_edges_lists();

        vector<vector<vector<vector<int>>>> route(n, vector<vector<vector<int>>>(n));
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                route[i][j] = yen_algorithm(edges_lists, i, j, k_max);
                route[j][i] = route[i][j];
                for (int k = 0; k < k_max; k++) {
                    reverse(route[j][i][k].begin(), route[j][i][k].end());
                }
            }
        }
        return route;
    }
};

struct RSAQuery {
    int from = -1;
    int to = -1;
    int width = 0;

    RSAQuery() {}
    RSAQuery(int _from, int _to, int _width) : from(_from), to(_to), width(_width) {}

    void output_query() const {
        cout << "(" << from << "->" << to << ", size:" << width << ")" << endl;
    }
};

struct RSAOnlineQuery : RSAQuery {
    double arrival = 0;
    double duration = 0;

    RSAOnlineQuery() {}
    RSAOnlineQuery(int _from, int _to, int _width) : RSAQuery(_from, _to, _width) {}
    RSAOnlineQuery(double _arrival, double _duration, int _from, int _to, int _width) : arrival(_arrival), duration(_duration), RSAQuery(_from, _to, _width) {}

    void output_query() const {
        cout << "arrival:" << arrival << ", duration:" << duration << ", ";
        RSAQuery::output_query();
    }
};

struct RSAAssignment {
    int fs = -1, fe = -1;
    vector<int> path;

    RSAAssignment() {}
    RSAAssignment(int _fs, int _fe, vector<int> _path) : fs(_fs), fe(_fe), path(_path) {}

    void output_assignment(const Graph &graph, const RSAQuery &query) const {
        if (fs == -1) {
            cout << "null" << endl;
            return;
        }

        auto edges = graph.get_edges();
        int now = query.from;
        cout << now;
        for (auto e_id : path) {
            now = edges[e_id].next(now);
            cout << "-[" << e_id << "]->" << now;
        }
        cout << endl;
    }
};

struct RSAQueryGenerator {
    int vertices_num = 1;
    int slots_num = 1;
    double time_max = 1.0;
    double width_avg = 1.0, width_var = 1.0;
    double query_num_avg = 1.0;
    double query_time_avg = 1.0, query_time_var = 1.0;

    random_device rnd;

    RSAQueryGenerator() {}
    RSAQueryGenerator(int _vertices_num, int _slots_num, double _time_max, double _width_avg, double _width_var, double _query_num_avg, double _query_time_avg, double _query_time_var) :
        vertices_num(_vertices_num), slots_num(_slots_num), time_max(_time_max), width_avg(_width_avg), width_var(_width_var),
        query_num_avg(_query_num_avg), query_time_avg(_query_time_avg), query_time_var(_query_time_var) {}

    vector<RSAOnlineQuery> generate_random_queries() {
        normal_distribution<> w_dist(width_avg, width_var);
        normal_distribution<> duration_dist(query_time_avg, query_time_var);
        vector<RSAOnlineQuery> queries;

        double now_time = 0;
        while (true) {
            double p = (double)rnd() / ((double)UINT_MAX + 1);
            now_time += -log(1 - p) / query_num_avg;
            if (now_time > time_max) break;

            double q = (double)rnd() / ((double)UINT_MAX + 1);
            double duration = max(duration_dist(rnd), 0.0);
            int s = rnd() % vertices_num, t = rnd() % vertices_num;
            while (s == t) {
                s = rnd() % vertices_num, t = rnd() % vertices_num;
            }
            int w = min(max((int)round(w_dist(rnd)), 1), slots_num);
            queries.push_back(RSAOnlineQuery(now_time, duration, s, t, w));
        }
        return queries;
    }
};

struct RSABenchmarker {
protected:
    Graph graph;
    int slots_num = 0;
    vector<vector<int>> slots;
    vector<RSAOnlineQuery> queries;
    vector<RSAAssignment> solution;
    int all_queries_num = 0;
    int success_num = 0;
    vector<vector<int>> slots_usages;

public:
    RSABenchmarker() {}
    RSABenchmarker(Graph _graph, int _slots_num) : graph(_graph), slots_num(_slots_num) {}
    RSABenchmarker(Graph _graph, int _slots_num, vector<RSAOnlineQuery> _queries) : graph(_graph), slots_num(_slots_num), queries(_queries) {}

    void set_graph(Graph _graph) {
        graph = _graph;
    }

    Graph get_graph() const {
        return graph;
    }

    void set_slots_num(int _slots_num) {
        slots_num = _slots_num;
    }

    int get_slots_num() const {
        return slots_num;
    }

    vector<vector<int>> get_slots() const {
        return slots;
    }

    void set_queries(vector<RSAOnlineQuery> _queries) {
        queries = _queries;
    }

    void add_query(RSAOnlineQuery query) {
        queries.push_back(query);
    }

    vector<RSAOnlineQuery> get_queries() const {
        return queries;
    }

    vector<RSAAssignment> get_solution() const {
        return solution;
    }

    void output_solution() const {
        int sz = queries.size();
        if (solution.size() != sz) {
            cout << "solution size error" << endl;
        }
        for (int i = 0; i < sz; i++) {
            cout << "query_id: " << i << ", ";
            queries[i].output_query();
            solution[i].output_assignment(graph, queries[i]);
            cout << endl;
        }
    }

    int get_all_queries_num() const {
        return all_queries_num;
    }

    int get_success_num() const {
        return success_num;
    }

    double get_blocking_rate() const {
        if (all_queries_num == 0) return 1;
        return (double)(all_queries_num - success_num) / all_queries_num;
    }

    void initiate() {
        const int n = graph.get_vertices_num();
        const int m = graph.get_edges().size();

        all_queries_num = 0, success_num = 0;
        slots_usages = vector<vector<int>>(m, vector<int>(slots_num, false));
    }

    virtual void run() = 0;
};

struct RSASolver {
protected:
    Graph graph;
    int slots_num = 0;
    vector<vector<int>> slots;
    vector<RSAQuery> queries;
    vector<RSAAssignment> solution;

    random_device rnd;

public:
    RSASolver() {}
    RSASolver(Graph _graph) : graph(_graph), slots_num(0) {}
    RSASolver(Graph _graph, int _slots_num) : graph(_graph), slots_num(_slots_num) {}
    RSASolver(Graph _graph, int _slots_num, vector<RSAQuery> _queries) : graph(_graph), slots_num(_slots_num), queries(_queries) {}

    void set_graph(Graph _graph) {
        graph = _graph;
    }

    Graph get_graph() const {
        return graph;
    }

    void set_slots_num(int _slots_num) {
        slots_num = _slots_num;
    }

    int get_slots_num() const {
        return slots_num;
    }

    void set_slots(vector<vector<int>> _slots) {
        slots = _slots;
    }

    vector<vector<int>> get_slots() const {
        return slots;
    }

    void set_queries(vector<RSAQuery> _queries) {
        queries = _queries;
    }

    void add_query(RSAQuery query) {
        queries.push_back(query);
    }

    vector<RSAQuery> get_queries() const {
        return queries;
    }

    bool validate() const {
        int n = graph.get_vertices_num();
        auto edges = graph.get_edges();
        int m = edges.size();
        int q = queries.size();

        if (!(queries.size() == solution.size())) {
            cout << "error" << endl;
            return false;
        }

        vector<vector<int>> slots_usages(m, vector<int>(slots_num, false));
        for (int i = 0; i < q; i++) {
            auto query = queries[i];
            auto assignment = solution[i];

            if (assignment.fs == -1) continue;
            int now = query.from;
            for (auto e_id : assignment.path) {
                if (edges[e_id].from == now) {
                    now = edges[e_id].to;
                } else if (edges[e_id].to == now) {
                    now = edges[e_id].from;
                } else {
                    cout << "path error" << endl;
                    return false;
                }
            }
            if (!(now == query.to)) {
                cout << "path error" << endl;
                return false;
            }

            if (!(0 <= assignment.fs && assignment.fe - assignment.fs == query.width && assignment.fe <= slots_num)) {
                cout << "slot error" << endl;
                return false;
            }

            for (int e_id : assignment.path) {
                for (int l = assignment.fs; l < assignment.fe; l++) {
                    if (slots_usages[e_id][l]) {
                        cout << "kaburi" << endl;
                        return false;
                    }
                    slots_usages[e_id][l] = true;
                }
            }
        }

        cout << "correct solution" << endl;
        return true;
    }

    void output_solution() const {
        int n = graph.get_vertices_num();
        auto edges = graph.get_edges();
        int m = edges.size();
        int q = queries.size();

        vector<vector<int>> slots_usages = slots;
        for (int i = 0; i < q; i++) {
            cout << "[" << queries[i].from << " -> " << queries[i].to << ", width:" << queries[i].width << "]" << endl;
            cout << "length:" << solution[i].path.size() << ", slot:[" << solution[i].fs << ", " << solution[i].fe << ")" << endl;

            int now = queries[i].from;
            cout << now;
            for (auto e_id : solution[i].path) {
                now = edges[e_id].from + edges[e_id].to - now;
                cout << " -> " << now;
                for (int l = solution[i].fs; l < solution[i].fe; l++) {
                    slots_usages[e_id][l] = true;
                }
            }
            cout << endl;
        }
        for (auto v : slots_usages) {
            for (auto x : v) {
                cout << (x ? 'o' : 'x');
            }
            cout << endl;
        }
    }

    vector<RSAAssignment> get_solution() const {
        return solution;
    }

    vector<vector<int>> get_slots_usages() const {
        int n = graph.get_vertices_num();
        auto edges = graph.get_edges();
        int m = edges.size();
        int q = queries.size();

        vector<vector<int>> slots_usages = slots;
        for (int i = 0; i < q; i++) {
            int now = queries[i].from;
            for (auto e_id : solution[i].path) {
                now = edges[e_id].from + edges[e_id].to - now;
                for (int l = solution[i].fs; l < solution[i].fe; l++) {
                    slots_usages[e_id][l] = true;
                }
            }
        }

        return slots_usages;
    }
};

struct RSASolverKSPFF : RSASolver {
public:
    RSASolverKSPFF() : RSASolver() {}
    RSASolverKSPFF(Graph _graph) : RSASolver(_graph) {}
    RSASolverKSPFF(Graph _graph, int _slots_num) : RSASolver(_graph, _slots_num) {}
    RSASolverKSPFF(Graph _graph, int _slots_num, vector<RSAQuery> _queries) : RSASolver(_graph, _slots_num, _queries) {}

    void solve(const int k_max) {
        if (queries.size() != 1) {
            cout << "SolverKSPFF queries size must be 1" << endl;
        }
        const int from = queries.front().from;
        const int to = queries.front().to;
        const int width = queries.front().width;
        const auto edges_lists = graph.make_edges_lists();
        const auto route = graph.yen_algorithm(edges_lists, from, to, k_max);

        for (int k = 0; k < k_max; k++) {
            vector<int> available_slots(slots_num, 1);
            for (int e_id : route[k]) {
                const auto& sv = slots[e_id];
                for (int i = 0; i < slots_num; i++) {
                    if (sv[i]) available_slots[i] = 0;
                }
            }
            vector<int> sum(slots_num + 1, 0);
            for (int i = 1; i <= slots_num; i++) {
                sum[i] = sum[i - 1] + available_slots[i - 1];
            }
            for (int i = 0; i <= slots_num - width; i++) {
                if (sum[i + width] - sum[i] == width) {
                    solution = { RSAAssignment(i, i + width, route[k]) };
                    return;
                }
            }
        }

        solution = { RSAAssignment() };
    }
};

struct RSABenchmarkerKSPFF : RSABenchmarker {
private:
    int k_max = 0;

public:
    RSABenchmarkerKSPFF() {}
    RSABenchmarkerKSPFF(Graph _graph, int _slots_num) : RSABenchmarker(_graph, _slots_num) {}
    RSABenchmarkerKSPFF(int _k_max, Graph _graph, int _slots_num, vector<RSAOnlineQuery> _queries) : k_max(_k_max), RSABenchmarker(_graph, _slots_num, _queries) {}

    void set_k_max(int _k_max) {
        k_max = _k_max;
    }

    void run() {
        const int n = graph.get_vertices_num();
        const int m = graph.get_edges().size();

        priority_queue<pair<double, int>> release_heap;

        int iteration = queries.size();

        cout << setfill('0');
        for (int t = 0; t < queries.size(); t++) {
            double time_end = queries[t].arrival;
            while (!release_heap.empty()) {
                auto p = release_heap.top();
                if (-p.first >= time_end) break;
                release_heap.pop();

                int id = p.second;
                auto que = queries[id];
                auto asg = solution[p.second];

                for (auto e_id : asg.path) {
                    for (int l = asg.fs; l < asg.fe; l++) {
                        slots_usages[e_id][l] = false;
                    }
                }
            }

            all_queries_num++;
            RSASolverKSPFF rsa(graph, slots_num);
            rsa.set_slots(slots_usages);
            rsa.set_queries({queries[t]});
            rsa.solve(k_max);
            if (rsa.validate()) {
                slots_usages = rsa.get_slots_usages();
                vector<RSAAssignment> sol = rsa.get_solution();
                for (auto asg : sol) {
                    if (asg.fs != -1) success_num++;
                }
                solution.insert(solution.end(), sol.begin(), sol.end());
            } else { // validation error
                vector<RSAAssignment> sol(1);
                solution.insert(solution.end(), sol.begin(), sol.end());
            }
            if (t % 10 == 9) {
                cout << "KSPFF [" << setw(5) << t << "/" << setw(5) << iteration << "] " << "blocking_rate:" << get_blocking_rate() << endl;
            }

            release_heap.push({ -time_end - queries[t].duration, t });
        }
    }
};

struct RSASolverCPLEX : RSASolver {
public:
    RSASolverCPLEX() : RSASolver() {}
    RSASolverCPLEX(Graph _graph) : RSASolver(_graph) {}
    RSASolverCPLEX(Graph _graph, int _slots_num) : RSASolver(_graph, _slots_num) {}
    RSASolverCPLEX(Graph _graph, int _slots_num, vector<RSAQuery> _queries) : RSASolver(_graph, _slots_num, _queries) {}

    void solve(int k_max, double solver_time_limit) {
        int n = graph.get_vertices_num();
        auto edges = graph.get_edges();
        int m = edges.size();
        int q = queries.size();

        assert(slots.size() == slots_num);
        for (int ed = 0; ed < m; ed++) {
            assert(slots[ed].size() == slots_num);
        }

        vector<vector<int>> e(n);
        for (int i = 0; i < m; i++) {
            e[edges[i].from].push_back(i);
            e[edges[i].to].push_back(i);
        }

        auto route = graph.calc_route(k_max);
        const auto edges_lists = graph.make_edges_lists();

        // 初期実行可能解用に KSP-FF を解く
        vector<RSAAssignment> kspff_solution;
        vector<int> kspff_solution_k;
        vector<vector<int>> kspff_slots = slots;

        for (auto query : queries) {
            const int from = query.from;
            const int to = query.to;
            const int width = query.width;
            bool found = false;
            for (int k = 0; k < k_max; k++) {
                vector<int> available_slots(slots_num, 1);
                for (int e_id : route[from][to][k]) {
                    const auto& sv = kspff_slots[e_id];
                    for (int i = 0; i < slots_num; i++) {
                        if (sv[i]) available_slots[i] = 0;
                    }
                }
                vector<int> sum(slots_num + 1, 0);
                for (int i = 1; i <= slots_num; i++) {
                    sum[i] = sum[i - 1] + available_slots[i - 1];
                }
                for (int i = 0; i <= slots_num - width; i++) {
                    if (sum[i + width] - sum[i] == width) {
                        kspff_solution.emplace_back(i, i + width, route[from][to][k]);
                        kspff_solution_k.push_back(k);
                        for (int e_id : route[from][to][k]) {
                            for (int j = i; j < i + width; j++) {
                                kspff_slots[e_id][j] = true;
                            }
                        }
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }

            if (found) continue;
            kspff_solution.emplace_back();
            kspff_solution_k.push_back(k_max);
        }

        // edge_usage[i][j][ed] := クエリ i の j 番目のパスが辺 e を使用しているか？
        vector<vector<vector<int>>> edge_usage(q, vector<vector<int>>(k_max, vector<int>(m, false)));

        for (int i = 0; i < q; i++) {
            for (int j = 0; j < k_max; j++) {
                for (int ed : route[queries[i].from][queries[i].to][j]) {
                    edge_usage[i][j][ed] = true;
                }
            }
        }

        IloEnv env;
        IloModel model(env);
        IloNumVarArray vars(env);

        // 数理計画問題の変数 fs[i] := リクエスト i における波長の始点
        for (int i = 0; i < q; i++) {
            vars.add(IloIntVar(env, 0, slots_num - queries[i].width));
        }
        auto get_fs = [&](int i) {
            return vars[i];
        };

        // 数理計画問題の変数 k[i] := リクエスト i における選択したパス
        for (int i = 0; i < q; i++) {
            vars.add(IloIntVar(env, 0, k_max));
        }
        const int k_offset = q;
        auto get_k = [&](int i) {
            return vars[i + k_offset];
        };

        // 数理計画問題の変数 x[i][j][l] := リクエスト i においてパス j を波長 l で使うときのみ 1
        for (int i = 0; i < q; i++) {
            for (int j = 0; j <= k_max; j++) {
                for (int l = 0; l < slots_num; l++) {
                    vars.add(IloIntVar(env, 0, 1));
                }
            }
        }
        const int x_offset = q + k_offset;
        auto get_x = [&](int i, int j, int l) {
            return vars[i * (k_max + 1) * slots_num + j * slots_num + l + x_offset];
        };

        // 初期解を設定
        IloNumVarArray start_vars(env);
        IloNumArray start_vals(env);
        for (int i = 0; i < q; i++) {
            auto asg = kspff_solution[i];
            int k = kspff_solution_k[i];
            int fs = (k == k_max ? 0 : asg.fs);
            start_vars.add(get_fs(i));
            start_vals.add(fs);
            start_vars.add(get_k(i));
            start_vals.add(k);
            for (int j = 0; j < k_max + 1; j++) {
                for (int l = 0; l < slots_num; l++) {
                    start_vars.add(get_x(i, j, l));
                    if (j == k && fs <= l && l < fs + queries[i].width) {
                        start_vals.add(1);
                    } else {
                        start_vals.add(0);
                    }
                }
            }
        }

        // 目的関数
        IloExpr min_exp(env);
        for (int i = 0; i < q; i++) {
            // 辺数 * 割当失敗数
            min_exp += m * get_x(i, k_max, 0);
        }
        for (int i = 0; i < q; i++) {
            for (int j = 0; j < k_max; j++) {
                int edge_num = route[queries[i].from][queries[i].to][j].size();
                int num = (slots_num + queries[i].width - 1) / queries[i].width;
                for (int l = 1; l < num; l++) {
                    // パスが含む辺の数 * スロットの使用コスト
                    min_exp += edge_num * (double)l / num * get_x(i, j, l * queries[i].width);
                }
            }
        }
        model.add(IloMinimize(env, min_exp));

        for (int i = 0; i < q; i++) {
            // 経路が存在しない j は使用しない
            for (int j = 0; j < k_max; j++) {
                if (route[queries[i].from][queries[i].to][j].size() == 0) {
                    for (int l = 0; l < slots_num; l++) {
                        model.add(get_x(i, j, l) == 0);
                    }
                }
            }

            // x[i][j][l] = 1 のとき、k[i] <= j <= k[i], fs[i] <= l < fs[i] + width
            for (int j = 0; j < k_max; j++) {
                for (int l = 0; l < slots_num; l++) {
                    model.add(get_k(i) - j <= (k_max + 1) * (1 - get_x(i, j, l)));
                    model.add(get_k(i) - j >= (k_max + 1) * (get_x(i, j, l) - 1));
                    model.add(get_fs(i) - l <= slots_num * (1 - get_x(i, j, l)));
                    model.add(get_fs(i) + queries[i].width - 1 - l >= slots_num * (get_x(i, j, l) - 1));
                }
            }

            // k[i] = k_max のとき、fs[i] = 0
            model.add(get_fs(i) <= slots_num * (k_max - get_k(i)));

            for (int l = 0; l < slots_num; l++) {
                if (l < queries[i].width) {
                    model.add(get_k(i) - k_max >= k_max * (get_x(i, k_max, l) - 1));
                } else {
                    model.add(get_x(i, k_max, l) == 0);
                }
            }

            // x[i] の和は query.width
            IloExpr exp_sum(env);
            for (int j = 0; j <= k_max; j++) {
                for (int l = 0; l < slots_num; l++) {
                    exp_sum += get_x(i, j, l);
                }
            }
            model.add(exp_sum == queries[i].width);
        }

        // リクエスト同士で同じ辺の同じ波長を使わない
        for (int ed = 0; ed < m; ed++) {
            for (int l = 0; l < slots_num; l++) {
                IloExpr exp(env);
                for (int i = 0; i < q; i++) {
                    for (int j = 0; j < k_max; j++) {
                        if (edge_usage[i][j][ed]) {
                            exp += get_x(i, j, l);
                        }
                    }
                }
                if (slots[ed][l]) {
                    model.add(exp == 0);
                } else {
                    model.add(exp <= 1);
                }
            }
        }

        IloCplex cplex(model);
        IloNumArray vals(env);

        cplex.addMIPStart(start_vars, start_vals);
        start_vals.end();
        start_vars.end();

        cplex.setOut(env.getNullStream());

        cplex.setParam(IloCplex::Param::TimeLimit, solver_time_limit);
        // cplex.setParam(IloCplex::Param::Emphasis::MIP, CPX_MIPEMPHASIS_FEASIBILITY);
        cplex.solve();
        cplex.getValues(vals, vars);

        // env.out() << "Solution Value = " << cplex.getObjValue() << endl;
        // env.out() << "Values = " << vals << endl;

        for (int i = 0; i < q; i++) {
            int fs = vals[i];
            int k = vals[i + k_offset];
            if (k != k_max) {
                solution.push_back(RSAAssignment(fs, fs + queries[i].width, route[queries[i].from][queries[i].to][k]));
            } else {
                solution.push_back(RSAAssignment(-1, -1, {}));
            }
        }

        env.end();
    }
};

struct RSABenchmarkerCPLEX : RSABenchmarker {
private:
    int iteration = 0;
    double time_max = 0.0;
    double solver_time_span = 1.0;
    double solver_time_limit = 0.0;
    int k_max = 0;

public:
    RSABenchmarkerCPLEX() {}
    RSABenchmarkerCPLEX(Graph _graph, int _slots_num) : RSABenchmarker(_graph, _slots_num) {}
    RSABenchmarkerCPLEX(Graph _graph, int _slots_num, vector<RSAOnlineQuery> _queries) : RSABenchmarker(_graph, _slots_num, _queries) {}
    RSABenchmarkerCPLEX(double _time_max, double _solver_time_span, double _solver_time_limit, int _k_max, Graph _graph, int _slots_num, vector<RSAOnlineQuery> _queries) :
        iteration(ceil(_time_max / _solver_time_span)), time_max(_time_max), solver_time_span(_solver_time_span), solver_time_limit(_solver_time_limit),
        k_max(_k_max), RSABenchmarker(_graph, _slots_num, _queries) {}

    void set_solver_properties(double _time_max, double _solver_time_span, double _solver_time_limit, int _k_max) {
        iteration = ceil(_time_max / _solver_time_span);
        time_max = _time_max;
        solver_time_span = _solver_time_span;
        solver_time_limit = _solver_time_limit;
        k_max = _k_max;
    }
    
    void run() {
        const int n = graph.get_vertices_num();
        const int m = graph.get_edges().size();

        priority_queue<pair<double, int>> release_heap;

        cout << setfill('0');
        for (int t = 0; t < iteration; t++) {
            double time_start = t * solver_time_span, time_end = (t + 1) * solver_time_span;
            while (!release_heap.empty()) {
                auto p = release_heap.top();
                if (-p.first >= time_end) break;
                release_heap.pop();

                int id = p.second;
                auto que = queries[id];
                auto asg = solution[p.second];

                for (auto e_id : asg.path) {
                    for (int l = asg.fs; l < asg.fe; l++) {
                        slots_usages[e_id][l] = false;
                    }
                }
            }

            vector<int> query_ids;
            vector<RSAQuery> qs;
            for (int i = 0; i < queries.size(); i++) {
                auto q = queries[i];
                if (time_start <= q.arrival && q.arrival < time_end) {
                    qs.push_back(q);
                    query_ids.push_back(i);
                }
            }

            all_queries_num += qs.size();
            RSASolverCPLEX rsa(graph, slots_num);
            rsa.set_slots(slots_usages);
            rsa.set_queries(qs);
            rsa.solve(k_max, solver_time_limit);
            if (rsa.validate()) {
                slots_usages = rsa.get_slots_usages();
                vector<RSAAssignment> sol = rsa.get_solution();
                for (auto asg : sol) {
                    if (asg.fs != -1) success_num++;
                }
                solution.insert(solution.end(), sol.begin(), sol.end());
            } else { // validation error
                vector<RSAAssignment> sol(qs.size());
                solution.insert(solution.end(), sol.begin(), sol.end());
            }
            cout << "CPLEX [" << setw(5) << t << "/" << setw(5) << iteration << "] " << "blocking_rate:" << get_blocking_rate() << endl;

            for (auto id : query_ids) {
                release_heap.push({ -time_end - queries[id].duration, id });
            }
        }
    }
};

void output_route(Graph& graph, vector<vector<vector<vector<int>>>>& route, int k_max) {
    const int n = graph.get_vertices_num();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            cout << "[" << i << "->" << j << "]" << endl;
            for (int k = 0; k < k_max; k++) {
                cout << k << ":";
                if (route[i][j][k].size() == 0) {
                    cout << "null" << endl;
                    continue;
                }
                int now = i;
                for (auto e : route[i][j][k]) {
                    cout << " " << now;
                    now = graph.get_edges()[e].next(now);
                }
                cout << " " << now << endl;
            }
        }
    }
}

void ksp_test() {
    Graph graph(4, {
        Edge(0, 1, 2),
        Edge(0, 2, 3),
        Edge(1, 2, 2),
        Edge(1, 3, 3),
        Edge(2, 3, 4),
    });
    const int k_max = 10;
    auto route = graph.calc_route(k_max);

    output_route(graph, route, k_max);
}

int main() {
    cout << std::fixed << std::setprecision(10) << endl;

    Graph nsf_graph(14, {
        Edge(0, 1, 1050),
        Edge(0, 2, 1500),
        Edge(0, 7, 2400),
        Edge(1, 2, 600),
        Edge(1, 3, 750),
        Edge(2, 5, 1800),
        Edge(3, 4, 600),
        Edge(3, 10, 1950),
        Edge(4, 5, 1200),
        Edge(4, 6, 600),
        Edge(5, 9, 1050),
        Edge(5, 13, 1800),
        Edge(6, 7, 750),
        Edge(6, 9, 1350),
        Edge(7, 8, 750),
        Edge(8, 9, 750),
        Edge(8, 11, 300),
        Edge(8, 12, 300),
        Edge(10, 11, 600),
        Edge(10, 12, 750),
        Edge(11, 13, 300),
        Edge(12, 13, 150),
    });
    const int n = nsf_graph.get_vertices_num();
    const int m = nsf_graph.get_edges().size();

    const double time_max = 1000;
    const int slots_num = 50, k_max = 6;
    const double width_avg = 5.0, width_var = 2.0;
    const double query_num_avg = 2.0;
    const double query_time_avg = 20.0, query_time_var = 10.0;
    const double solver_time_span = 5.0;
    const double solver_time_limit = 5.0;

    RSAQueryGenerator query_generator(n, slots_num, time_max, width_avg, width_var, query_num_avg, query_time_avg, query_time_var);
    vector<RSAOnlineQuery> queries = query_generator.generate_random_queries();

    
    RSABenchmarkerCPLEX bench_cplex(time_max, solver_time_span, solver_time_limit, k_max, nsf_graph, slots_num, queries);
    bench_cplex.initiate();
    bench_cplex.run();
    
    bench_cplex.output_solution();

    RSABenchmarkerKSPFF bench_kspff(k_max, nsf_graph, slots_num, queries);
    bench_kspff.initiate();
    bench_kspff.run();

    bench_kspff.output_solution();

    cout << "all_queries_num:" << bench_cplex.get_all_queries_num() << ", success_num:" << bench_cplex.get_success_num() << endl;
    cout << "blocking_rate:" << bench_cplex.get_blocking_rate() << endl;
    cout << "all_queries_num:" << bench_kspff.get_all_queries_num() << ", success_num:" << bench_kspff.get_success_num() << endl;
    cout << "blocking_rate:" << bench_kspff.get_blocking_rate() << endl;

    return 0;
}