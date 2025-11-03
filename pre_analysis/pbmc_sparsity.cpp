#include <bits/stdc++.h>
using namespace std;

#define ll long long
#define vs vector<string>
#define vi vector<int>
#define vvi vector<vi>
#define pb push_back
#define P(n) (cout << #n << ": " << (n) << '\n')

struct scRNA_matrix {
    vvi mat;
    vs cells;
    vs genes;
    vs labels;
} dat;

vector<string> split_csv_line(const string& line) {
    vector<string> tokens;
    string item;
    stringstream ss(line);
    while (getline(ss, item, ',')) tokens.pb(item);
    return tokens;
}

template <class T>
vector<vector<T>> transpose(const vector<vector<T>>& A) {
    if (A.empty()) return {};
    size_t R = A.size(), C = A[0].size();
    vector<vector<T>> B(C, vector<T>(R));
    for (size_t i = 0; i < R; ++i)
        for (size_t j = 0; j < C; ++j)
            B[j][i] = A[i][j];
    return B;
}

void Read_pbmc(const string& file_matrix, const string& file_labels) {
    ifstream fin(file_matrix);
    if (!fin) {
        cerr << "Failed to open matrix file: " << file_matrix << "\n";
        exit(1);
    }
    string line;
    if (!getline(fin, line)) {
        cerr << "Matrix file is empty.\n";
        exit(1);
    }
    auto head = split_csv_line(line);
    if (head.size() < 2) {
        cerr << "Matrix header has <2 columns (expect: gene_name, cell1, cell2, ...).\n";
        exit(1);
    }
    dat.cells.assign(head.begin() + 1, head.end());
    int size_cells = (int)dat.cells.size();
    P(size_cells);

    vvi rows; rows.reserve(50000);
    vs genes; genes.reserve(50000);

    int line_no = 1;
    while (getline(fin, line)) {
        ++line_no;
        if (line.empty()) continue;
        auto toks = split_csv_line(line);
        if ((int)toks.size() != size_cells + 1) {
            cerr << "Bad line at " << line_no << ": expect " << (size_cells + 1)
                 << " columns, got " << toks.size() << ".\n";
            exit(1);
        }
        genes.pb(toks[0]);
        vi tmp(size_cells);
        for (int i = 0; i < size_cells; ++i) {
            if (toks[i + 1].empty()) tmp[i] = 0;
            else {
                
                tmp[i] = stoi(toks[i + 1]);
            }
        }
        rows.pb(move(tmp));
    }
    fin.close();

    dat.genes = move(genes);
    int size_genes = (int)rows.size();
    P(size_genes);

    
    dat.mat = transpose(rows);

    ifstream fin2(file_labels);
    if (!fin2) {
        cerr << "Warning: failed to open label file: " << file_labels << " (labels will be empty)\n";
        dat.labels.assign(dat.cells.size(), "");
        return;
    }
    string l2;
    if (!getline(fin2, l2)) {
        cerr << "Warning: empty label file.\n";
        dat.labels.assign(dat.cells.size(), "");
        return;
    }
    unordered_map<string, string> label_map;
    while (getline(fin2, l2)) {
        if (l2.empty()) continue;
        auto tt = split_csv_line(l2);
        if (tt.empty()) continue;
        string cell = tt[0];
        string lab  = (tt.size() >= 2 ? tt[1] : "");
        label_map.emplace(move(cell), move(lab));
    }
    fin2.close();

    dat.labels.resize(dat.cells.size());
    int hit = 0;
    for (size_t i = 0; i < dat.cells.size(); ++i) {
        auto it = label_map.find(dat.cells[i]);
        if (it != label_map.end()) {
            dat.labels[i] = it->second;
            ++hit;
        } else {
            dat.labels[i] = "";
        }
    }
    cout << "Label coverage: " << hit << " / " << dat.cells.size()
         << " = " << fixed << setprecision(4) << (double)hit / max<size_t>(1, dat.cells.size()) << "\n";
}

double quantile_from_sorted(const vector<int>& v, double q) {
    if (v.empty()) return 0.0;
    if (q <= 0) return v.front();
    if (q >= 1) return v.back();
    double idx = q * (v.size() - 1);
    size_t lo = (size_t)floor(idx), hi = (size_t)ceil(idx);
    if (lo == hi) return v[lo];
    double w = idx - lo;
    return v[lo] * (1.0 - w) + v[hi] * w;
}


void ReportSparsity(const vvi& A, int thr = 0) {
    if (A.empty() || A[0].empty()) {
        cout << "Empty matrix.\n";
        return;
    }
    const int R = (int)A.size();       // cells
    const int C = (int)A[0].size();    // genes
    for (int i = 1; i < R; ++i) {
        if ((int)A[i].size() != C) {
            cerr << "Jagged matrix at row " << i << ".\n";
            exit(1);
        }
    }

    long long nnz = 0;
    vector<int> row_nnz(R, 0), col_nnz(C, 0);

    for (int i = 0; i < R; ++i) {
        const auto& row = A[i];
        for (int j = 0; j < C; ++j) {
            if (row[j] > thr) {
                ++nnz;
                ++row_nnz[i];
                ++col_nnz[j];
            }
        }
    }

    const long long total = 1LL * R * C;
    const double density  = (double)nnz / (double)total;
    const double sparsity = 1.0 - density;
    const double avg_row_nnz = (double)nnz / (double)R;
    const double avg_col_nnz = (double)nnz / (double)C;

    auto row_sorted = row_nnz; sort(row_sorted.begin(), row_sorted.end());
    auto col_sorted = col_nnz; sort(col_sorted.begin(), col_sorted.end());

    cout << fixed << setprecision(6);
    cout << "Cells (R): " << R << ", Genes (C): " << C << '\n';
    cout << "Total entries: " << total << ", NNZ (> " << thr << "): " << nnz << '\n';
    cout << "Density (NNZ/total): " << density << '\n';
    cout << "Sparsity (1 - density): " << sparsity << '\n';

    auto rmin = row_sorted.front();
    auto rmed = quantile_from_sorted(row_sorted, 0.5);
    auto rp90 = quantile_from_sorted(row_sorted, 0.9);
    auto rmax = row_sorted.back();

    auto cmin = col_sorted.front();
    auto cmed = quantile_from_sorted(col_sorted, 0.5);
    auto cp90 = quantile_from_sorted(col_sorted, 0.9);
    auto cmax = col_sorted.back();

    cout << "Avg NNZ per cell (s2_row): " << avg_row_nnz
         << "   [min/med/p90/max: " << rmin << '/' << rmed << '/' << rp90 << '/' << rmax << "]\n";
    cout << "Avg NNZ per gene (s2_col): " << avg_col_nnz
         << "   [min/med/p90/max: " << cmin << '/' << cmed << '/' << cp90 << '/' << cmax << "]\n";
}

int main(int argc, const char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string file_matrix = "./datasets/PBMC-Zheng2017/PBMC_SC1.csv";
    string file_labels = "./datasets/PBMC-Zheng2017/PBMCLabels_SC1ClusterLabels.csv";
    int thr = 0;

    if (argc >= 2) file_matrix = argv[1];
    if (argc >= 3) file_labels = argv[2];
    if (argc >= 4) thr         = stoi(argv[3]);

    cout << "Matrix: " << file_matrix << "\n";
    cout << "Labels: " << file_labels << "\n";
    cout << "Nonzero threshold: value > " << thr << " treated as nonzero.\n\n";

    Read_pbmc(file_matrix, file_labels);
    ReportSparsity(dat.mat, thr);

    return 0;
}
