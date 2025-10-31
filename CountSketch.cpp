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
#define vvpii vector<vector<pii>>
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
vvi Countsketch_cols(const vvi& mat, uint32_t s);
vvpii BuildSNNGraph(const vvi& points, int k);
vi LouvainClustering(const vvpii& snnGraph);


int main(int argc, const char* argv[]) {
    string file_matrix, file_labels;
    file_matrix = "./datasets/PBMC-Zheng2017/PBMC_SC1.csv";
    file_labels = "./datasets/PBMC-Zheng2017/PBMCLabels_SC1ClusterLabels.csv";
    // Read into global dat
    Read_pbmc(file_matrix, file_labels);
    
    // CountSketch
    int s = 20;
    auto dat_countsketch = Countsketch_cols(dat.mat, s);
    // for (auto& r : dat_countsketch) { for (auto x : r) cout << x << " "; cout << "\n"; }
    
    // Build KNN
    int k = 15;
    vvpii snnGraph = BuildSNNGraph(dat_countsketch, k);

    // Clustering
    vi result_CountSketch = LouvainClustering(snnGraph);
    
    return 0;
}

/*
    Split a line of tokens separate by "," to an array.
    Used for Gene names.
*/
vector<string> split_csv_line(const string& line) {
    vector<string> tokens;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ',')) tokens.pb(item);
    return tokens;
}

/*
    Transpose the matrix.
*/
template <class T>
vector<vector<T>> transpose(const vector<vector<T>>& A) {
    if (A.empty()) return {};
    size_t R = A.size(), C = A[0].size();
    vector<vector<T>> B(C, vector<T>(R));
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
    for (size_t i = 0; i < dat.cells.size(); ++i) {
        auto it = map_label.find(dat.cells[i]);
        dat.labels[i] = (it == map_label.end() ? "" : it->second);
    }
    return;
}

/*
    Column hash.
*/
static inline uint32_t h_mod(uint64_t x, uint32_t s){
    const uint64_t a = 11400714819323198485ull;
    const uint64_t b = 0x9e3779b97f4a7c15ull;
    return (uint32_t)((a * x + b) % s);
}

/*
    Sign hash.
*/
static inline int xi(uint64_t x){
    const uint64_t a = 0xbf58476d1ce4e5b9ull;
    const uint64_t b = 0x94d049bb133111ebull;
    return ((a * x + b) >> 63) ? +1 : -1;
}

/*
    CountSketch.
    n*m -> n*s.
*/
vvi Countsketch_cols(const vvi& mat, uint32_t s) {
    int n = mat.size(), m = mat[0].size();
    vvi out(n, vi (s, 0));
    for (int i = 0; i < n; ++ i) {
        const auto& row = mat[i];
        for (int j = 0; j < m; ++ j) {
            int v = row[j];
            if (!v) continue; // Sparsity
            int c = h_mod(j, s);
            int sgn = xi(j);
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
    for (int i = 0; i < a.size(); ++ i) {
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
            if (pq.size() > k)
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
            vector<int> intersection;
            set_intersection(
                kNNSet[i].begin(), kNNSet[i].end(),
                kNNSet[j].begin(), kNNSet[j].end(),
                back_inserter(intersection)
            );
            int snnWeight = intersection.size();
            if (snnWeight > 0) {
                snnGraph[i].pb({j, snnWeight});
                snnGraph[j].pb({i, snnWeight});    
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
    igraph_create(&g, &edges, 0, IGRAPH_UNDIRECTED);
    
    igraph_community_multilevel(&g, &weights, 1.0, &membership, NULL, NULL);   

    cout << "Dynamic Louvain Clustering Results:" << '\n';
    int num_vertices = igraph_vcount(&g);
    vector<int> memb(num_vertices);
    for (int i = 0; i < num_vertices; ++ i) {
        memb[i] = (int)VECTOR(membership)[i];
        cout << i << " in cluster " << VECTOR(membership)[i] << '\n';
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&g);
    return memb;
}