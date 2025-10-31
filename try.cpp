#include <bits/stdc++.h>
#define ll long long
#define vs vector <string>
#define vi vector <int>
#define vvi vector <vi>
#define pb push_back
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

vector<string> split_csv_line(const string& line);
void read_pbmc(string file_matrix, string file_labels);
vector<vector<ll>> countsketch_cols(const vvi& mat, uint32_t s);


int main(int argc, const char* argv[]) {
    string file_matrix, file_labels;
    file_matrix = "./datasets/PBMC-Zheng2017/PBMC_SC1.csv";
    file_labels = "./datasets/PBMC-Zheng2017/PBMCLabels_SC1ClusterLabels.csv";
    // Read into global dat
    read_pbmc(file_matrix, file_labels);
    
    int k = 20;
    auto AS = countsketch_cols(dat.mat, k);
    cout << "\nAS: " << AS.size() << "x" << (AS.empty()?0:AS[0].size()) << "\n";
    for (auto& r : AS){ for (auto x : r) cout << x << " "; cout << "\n"; }        

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
void read_pbmc(string file_matrix, string file_labels) {
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

static inline uint32_t next_pow2(uint32_t x){
    if (x <= 1) return 1u;
    --x; x |= x>>1; x |= x>>2; x |= x>>4; x |= x>>8; x |= x>>16;
    return x+1;
}
static inline uint32_t h(uint64_t x, uint32_t t){
    const uint64_t a = 11400714819323198485ull;
    const uint64_t b = 0x9e3779b97f4a7c15ull;
    return (uint32_t)((a * x + b) >> (64 - t));
}
static inline int xi(uint64_t x){
    const uint64_t a = 0xbf58476d1ce4e5b9ull;
    const uint64_t b = 0x94d049bb133111ebull;
    return ((a * x + b) >> 63) ? +1 : -1;
}

static inline uint32_t ilog2_pow2(uint32_t S){
    uint32_t t = 0;
    while ((1u << t) < S)
        ++ t;
    return t;
}

vector<vector<ll>> countsketch_cols(const vvi& mat, uint32_t s) {
    const uint32_t m = (uint32_t) mat.size();
    const uint32_t n = m ? (uint32_t) mat[0].size() : 0;

    uint32_t S = next_pow2(max(1u, s));
    uint32_t t = ilog2_pow2(S);

    vector<vector<ll>> out(m, vector<ll>(S, 0));

    for (uint32_t i = 0; i < m; ++i){
        const auto& row = mat[i];
        for (uint32_t j = 0; j < n; ++j){
            int v = row[j];
            if (v == 0) continue;
            uint32_t c = h(j, t);
            int sgn = xi(j);
            out[i][c] += (ll) sgn * (ll) v;
        }
    }
    return out;
}