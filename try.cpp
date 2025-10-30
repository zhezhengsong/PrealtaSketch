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

// vector<string> split_csv_line(const string& line);
void read_pbmc(string file_matrix, string file_labels);


int main(int argc, const char* argv[]) {
    string file_matrix, file_labels;
    file_matrix = "./datasets/PBMC-Zheng2017/PBMC_SC1.csv";
    file_labels = "./datasets/PBMC-Zheng2017/PBMCLabels_SC1ClusterLabels.csv";
    read_pbmc(file_matrix, file_labels);
    

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
    // for (auto x : dat.mat[22])
    //     cout << x << ' ';
    // cout << '\n';
    // P(dat.mat.size());
    dat.mat = transpose(dat.mat);
    // P(dat.mat.size());
    // for (auto x : dat.mat[22])
    //     cout << x << ' ';
    // cout << '\n';
    
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
    // for (int i = 0; i < dat.cells.size(); ++ i) {
    //     cout << dat.cells[i] << ' ' << dat.labels[i] << '\n';
    // }
    // cout << "\n";
    return;
}