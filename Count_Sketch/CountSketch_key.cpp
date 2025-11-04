// Full CountSketch.
// Multiple key version.
#include <bits/stdc++.h>
extern "C" {
  #include <igraph/igraph.h>
}
#define ll long long
#define pii pair <int, int>
#define pli pair <ll, int>
#define vs vector <string>
#define vi vector <int>
#define vvi vector <vi>
#define vvpii vector <vector <pii>>
#define pb push_back
#define fi first
#define se second
#define P(n) (cout << #n << ": " << (n) << '\n')
using namespace std;

// rows = genes, cols = cells
// labels aligned to cells
struct scRNA_matrix {
    vvi mat;
    vs cells;
    vs genes;
    vs labels;
} dat;

void Read_pbmc(string file_matrix, string file_labels);
uint64_t make_random_key();
vvi Countsketch_cols(const vvi& mat, uint32_t s, uint64_t key);
vvpii BuildSNNGraph(const vvi& points, int k);
vi LouvainClustering(const vvpii& snnGraph);
double Clustering_Accuracy(const vi& y_true, const vi& y_pred);

int main(int argc, const char* argv[]) {
    string file_matrix, file_labels;
    file_matrix = "./datasets/PBMC-Zheng2017/PBMC_SC1.csv";
    file_labels = "./datasets/PBMC-Zheng2017/PBMCLabels_SC1ClusterLabels.csv";
    // Read into global dat
    Read_pbmc(file_matrix, file_labels);

    uint64_t key = make_random_key();
    cout << "CountSketch key = 0x" << hex << key << dec << "\n";
    
    // CountSketch
    int s = 1000;
    auto dat_countsketch = Countsketch_cols(dat.mat, s, key);
    // for (auto& r : dat_countsketch) { for (auto x : r) cout << x << " "; cout << "\n"; }
    
    // Build KNN
    int k = 20;
    vvpii snnGraph = BuildSNNGraph(dat_countsketch, k);

    // Clustering
    vi result_CountSketch = LouvainClustering(snnGraph);

    // Calculate the accuracy
    vi y_true;
    for (auto x : dat.labels)
        y_true.pb(stof(x));
    // cout << "y_true:\n";
    // for (auto x : y_true) {
    //     cout << x << ' ';
    // }
    // cout << '\n';
    double acc = Clustering_Accuracy(y_true, result_CountSketch);
    cout << fixed << setprecision(6) << "Clustering Accuracy = " << acc << "\n";
    return 0;
}

/*
    Split a line of tokens separate by "," to an array.
    Used for Gene names.
*/
vs split_csv_line(const string& line) {
    vs tokens;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ','))
        tokens.pb(item);
    return tokens;
}

/*
    Transpose the matrix.
*/
template <class T>
vector <vector <T>> transpose(const vector <vector <T>>& A) {
    if (A.empty()) return {};
    size_t R = A.size(), C = A[0].size();
    vector <vector <T>> B(C, vector <T>(R));
    for (size_t i = 0; i < R; ++ i)
        for (size_t j = 0; j < C; ++ j)
            B[j][i] = A[i][j];
    return B;
}

/*
    Read from pbmc into dat.
    This data is from "the preprint".
*/
void Read_pbmc(string file_matrix, string file_labels) {
    // Read matrix
    ifstream fin(file_matrix);
    string line;
    getline(fin, line);
    auto tmp_cells = split_csv_line(line);
    dat.cells.assign(tmp_cells.begin() + 1, tmp_cells.end());
    int size_cells = dat.cells.size();
    P(size_cells);
    
    vvi rows;
    vs genes;
    while (getline(fin, line)) {
        auto tmp_row = split_csv_line(line);
        genes.pb(tmp_row[0]);
        vi tmp_v(size_cells);
        for (int i = 0; i < size_cells; ++ i)
            tmp_v[i] = stof(tmp_row[i + 1]);
        dat.mat.pb(tmp_v); // move
    }
    int size_genes = dat.mat.size();
    dat.genes = genes;
    P(size_genes);
    dat.mat = transpose(dat.mat);

    // Read labels
    ifstream fin2(file_labels);
    unordered_map <string, string> map_label;
    string l2;
    getline(fin2, l2);
    while (getline(fin2, l2)) {
        auto tt = split_csv_line(l2);
        map_label.emplace(tt[0], tt[1]);
    }
    dat.labels.resize(dat.cells.size());
    for (size_t i = 0; i < dat.cells.size(); ++ i) {
        auto it = map_label.find(dat.cells[i]);
        dat.labels[i] = (it == map_label.end() ? "" : it->second);
    }
    return;
}

uint64_t splitmix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/*
    Generate a random key.
*/
uint64_t make_random_key() {
    random_device rd;
    uint64_t a = (uint64_t(rd()) << 32) ^ uint64_t(rd());
    uint64_t b = (uint64_t(rd()) << 32) ^ uint64_t(rd());
    uint64_t t = (uint64_t)chrono::high_resolution_clock::now().time_since_epoch().count();
    return splitmix64(a ^ (b + 0x9e3779b97f4a7c15ULL) ^ t);
}

/*
    Column hash.
*/
uint32_t h_mod(uint64_t x, uint32_t s, uint64_t k1) {
    uint64_t y = splitmix64(x ^ k1);
    return (uint32_t)(y % s);
}

/*
    Sign hash.
*/
int xi(uint64_t x, uint64_t k2) {
    uint64_t y = splitmix64(x ^ k2);
    return (y >> 63) ? +1 : -1;
}

/*
    CountSketch.
    n*m -> n*s.
*/
vvi Countsketch_cols(const vvi& mat, uint32_t s, uint64_t key) {
    int n = mat.size(), m = mat[0].size();
    vvi out(n, vi (s, 0));
    uint64_t k1 = splitmix64(key);
    uint64_t k2 = splitmix64(key ^ 0x7f4a7c159e3779b9ULL);
    for (int i = 0; i < n; ++ i) {
        const auto& row = mat[i];
        for (int j = 0; j < m; ++ j) {
            int v = row[j];
            if (!v) continue; // Sparsity
            int c = h_mod((uint64_t)j, s, k1);
            int sgn = xi((uint64_t)j, k2);
            out[i][c] += sgn * v;
        }
    }
    return out;
}

/*
    Calculate the euclidean distance.
*/
ll euclideanDistance(const vi& a, const vi& b) { // int?
    ll sum = 0;
    for (size_t i = 0; i < a.size(); ++ i) {
        ll tmp = (a[i] - b[i]) * (a[i] - b[i]);
        sum += tmp;
    }
    return sum;
}

/*
    Build SNN from a vvi.
*/
vvpii BuildSNNGraph(const vvi& points, int k) {
    int N = points.size();
    
    // KNN
    vvi kNNGraph(N); // kNNGraph[i] is the KNN of node i.
    vector <set <int>> kNNSet(N); // Set ver of kNNGraph.
    for (int i = 0; i < N; ++ i) {
        priority_queue<pli> pq;
        for (int j = 0; j < N; ++ j) {
            if (i == j) continue;
            ll dist = euclideanDistance(points[i], points[j]);
            pq.push({dist, j});
            if (pq.size() > (size_t) k)
                pq.pop();
        }
        // cout << "point " << i << "'s k-NN: ";
        while (!pq.empty()) {
            int neighborIndex = pq.top().se;
            kNNGraph[i].push_back(neighborIndex);
            kNNSet[i].insert(neighborIndex);
            // cout << neighborIndex << " ";
            pq.pop();
        }
        // cout << '\n';
    }

    // SNN
    vvpii snnGraph(N);
    for (int i = 0; i < N; ++ i) {
        for (int j = i + 1; j < N; ++ j) {
            vi intersection;
            set_intersection(
                kNNSet[i].begin(), kNNSet[i].end(),
                kNNSet[j].begin(), kNNSet[j].end(),
                back_inserter(intersection)
            );
            int snnWeight = intersection.size();
            if (snnWeight > 0) {
                snnGraph[i].pb({j, snnWeight});
                // snnGraph[j].pb({i, snnWeight}); // Important fix?
                // cout << " find edge:(" << i << ", " << j << "), SNN weight = " << snnWeight << '\n';
            }
        }
    }
    return snnGraph;
}

/*
    My add_edge for igraph.
*/
void My_add_edge(igraph_vector_int_t* edges, igraph_vector_t* weights, int from, int to, double weight) {
    igraph_vector_int_push_back(edges, from);
    igraph_vector_int_push_back(edges, to);
    igraph_vector_push_back(weights, weight);
}

/*
    Clustering on SNN using igraph.
*/
vi LouvainClustering(const vvpii& snnGraph) {
    igraph_set_error_handler(igraph_error_handler_abort);
    igraph_t g;
    igraph_vector_int_t edges;
    igraph_vector_t weights;  
    igraph_vector_int_t membership;
    igraph_vector_int_init(&edges, 0);
    igraph_vector_init(&weights, 0);
    igraph_vector_int_init(&membership, 0);

    int N = snnGraph.size();
    for (int i = 0; i < N; ++ i)
        for (auto x : snnGraph[i])
            My_add_edge(&edges, &weights, i, x.fi, x.se);
    igraph_create(&g, &edges, N, IGRAPH_UNDIRECTED);
    
    igraph_community_multilevel(&g, &weights, 1.2, &membership, NULL, NULL);   

    // cout << "Dynamic Louvain Clustering Results:" << '\n';
    int num_vertices = igraph_vcount(&g);
    vi memb(num_vertices);
    for (int i = 0; i < num_vertices; ++ i) {
        memb[i] = (int)VECTOR(membership)[i];
        // cout << i << " in cluster " << VECTOR(membership)[i] << '\n';
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&g);
    return memb;
}

/*
    Hungarian algo.
*/
vi hungarian_min_cost(const vvi& cost) {
    int n = cost.size();
    const int INF = numeric_limits<int>::max() / 4;

    vi u(n + 1, 0), v(n + 1, 0), p(n + 1, 0), way(n + 1, 0);
    for (int i = 1; i <= n; ++ i) {
        p[0] = i;
        int j0 = 0;
        vi minv(n + 1, INF);
        vector<char> used(n + 1, false);

        do {
            used[j0] = true;
            int i0 = p[j0], j1 = 0;
            int delta = INF;
            for (int j = 1; j <= n; ++ j)
                if (!used[j]) {
                    int cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
            }
            for (int j = 0; j <= n; ++ j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                }
                else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);

        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    vi assignment(n, -1);
    for (int j = 1; j <= n; ++ j)
        if (p[j] != 0) {
            assignment[p[j] - 1] = j - 1;
    }
    return assignment;
}

/*
    Calculate the clustering accuracy.
*/
double Clustering_Accuracy(const vi& y_true, const vi& y_pred) {
    if (y_true.size() != y_pred.size() || y_true.empty())
        return 0.0;

    unordered_map <int, int> map_true, map_pred;
    int rt = 0, rp = 0;
    for (int t : y_true)
        if (!map_true.count(t))
            map_true[t] = rt ++;
    for (int p : y_pred)
        if (!map_pred.count(p))
            map_pred[p] = rp ++;

    int k = max(rt, rp);
    if (k == 0) return 0.0;

    vvi C(k, vi (k, 0));
    for (size_t i = 0; i < y_true.size(); ++ i) {
        int r = map_true[y_true[i]];
        int c = map_pred[y_pred[i]];
        C[r][c] += 1;
    }

    int M = 0;
    for (int i = 0; i < k; ++ i)
        for (int j = 0; j < k; ++ j)
            M = max(M, C[i][j]);

    vvi cost(k, vi (k, 0));
    for (int i = 0; i < k; ++ i)
        for (int j = 0; j < k; ++ j)
            cost[i][j] = M - C[i][j];

    vi assign = hungarian_min_cost(cost);

    ll matched = 0;
    for (int i = 0; i < k; ++ i) {
        int j = assign[i];
        if (j >= 0)
            matched += C[i][j];
    }

    return (double)matched / (double)y_true.size();
}