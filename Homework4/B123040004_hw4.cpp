/*
 * HW4: Solving TSP with Ant Colony Optimization (ACO)
 * Final Version: Added "mean distance" output to match requirements.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <climits>
#include <algorithm>
#include <iomanip>
#include <random>

using namespace std;

// --- 全域隨機數產生器 ---
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0.0, 1.0);

// --- Point 類別 ---
class Point {
public:
    int PointNumber;
    int X;
    int Y;
    
    Point(int pointNumber = 0, int x = 0, int y = 0) {
        PointNumber = pointNumber;
        X = x;
        Y = y;
    }
};

struct Ant {
    vector<int> path;       // 路徑
    vector<bool> visited;   // 訪問標記
    double tour_length;     // 總長度

    Ant(int n) {
        visited.resize(n, false);
        tour_length = 0.0;
        path.reserve(n + 1);
    }

    void clear(int n) {
        path.clear();
        fill(visited.begin(), visited.end(), false);
        tour_length = 0.0;
    }
};

class ACO {
private:
    int num_ants;
    double alpha;
    double beta;
    double rho;
    double Q;
    int max_runs;
    
    vector<Point> cities;
    int n;
    
    vector<vector<double>> dist_matrix;
    vector<vector<double>> pheromone;
    vector<vector<double>> heuristic_cache;

    // 最佳解紀錄
    vector<int> global_best_path;
    double global_best_dist;
    
    // 新增：平均距離紀錄
    double mean_best_dist; 

public:
    ACO(const vector<Point>& pts) : cities(pts), n(pts.size()) {
        num_ants = 30;      
        alpha = 1.0;
        beta = 2.0;         
        rho = 0.1;
        Q = 100.0;
        max_runs = 30;      
        
        global_best_dist = 1e18;
        mean_best_dist = 0.0;

        // 計算距離矩陣
        dist_matrix.resize(n, vector<double>(n));
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                double dx = (double)cities[i].X - cities[j].X;
                double dy = (double)cities[i].Y - cities[j].Y;
                dist_matrix[i][j] = sqrt(dx*dx + dy*dy);
            }
        }

        // 預先計算啟發值
        heuristic_cache.resize(n, vector<double>(n));
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                if(i != j) {
                    double dist = dist_matrix[i][j];
                    if (dist < 1e-10) dist = 1e-10;
                    double eta = 1.0 / dist;
                    heuristic_cache[i][j] = pow(eta, beta);
                }
            }
        }
    }

    void init_pheromones() {
        pheromone.assign(n, vector<double>(n, 1.0));
    }

    void construct_solution(Ant& ant) {
        ant.clear(n);
        int start_node = rand() % n;
        ant.path.push_back(start_node);
        ant.visited[start_node] = true;

        int current = start_node;
        for(int step = 0; step < n - 1; ++step) {
            int next_node = select_next_node(ant, current);
            ant.path.push_back(next_node);
            ant.visited[next_node] = true;
            ant.tour_length += dist_matrix[current][next_node];
            current = next_node;
        }
        ant.tour_length += dist_matrix[current][start_node];
    }

    int select_next_node(const Ant& ant, int current) {
        double sum_prob = 0.0;
        for(int j=0; j<n; ++j) {
            if(!ant.visited[j]) {
                double tau = pheromone[current][j];
                double eta_pow_beta = heuristic_cache[current][j];
                double p = (alpha == 1.0) ? (tau * eta_pow_beta) : (pow(tau, alpha) * eta_pow_beta);
                sum_prob += p;
            }
        }

        double r = dis(gen) * sum_prob;
        double cum_prob = 0.0;
        for(int j=0; j<n; ++j) {
            if(!ant.visited[j]) {
                double tau = pheromone[current][j];
                double eta_pow_beta = heuristic_cache[current][j];
                double p = (alpha == 1.0) ? (tau * eta_pow_beta) : (pow(tau, alpha) * eta_pow_beta);
                cum_prob += p;
                if(r <= cum_prob) return j;
            }
        }
        for(int j=0; j<n; ++j) if(!ant.visited[j]) return j;
        return -1;
    }

    void update_pheromones(const vector<Ant>& ants) {
        double evaporation = 1.0 - rho;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                pheromone[i][j] *= evaporation;
            }
        }
        for(const auto& ant : ants) {
            double delta = Q / ant.tour_length;
            for(size_t i=0; i < ant.path.size() - 1; ++i) {
                int u = ant.path[i];
                int v = ant.path[i+1];
                pheromone[u][v] += delta;
                pheromone[v][u] += delta;
            }
            int u = ant.path.back();
            int v = ant.path[0];
            pheromone[u][v] += delta;
            pheromone[v][u] += delta;
        }
    }

    void run() {
        cout << "Dataset size: " << n << " cities." << endl;
        cout << "Max Evaluations per run: " << 10000 * n << endl;

        double sum_best_dists = 0;
        int max_evals = 10000 * n;
        vector<Ant> ants(num_ants, Ant(n));

        for(int run = 0; run < max_runs; ++run) {
            init_pheromones();
            double run_best_dist = 1e18;
            int eval_count = 0;

            cout << "Running " << setw(2) << run + 1 << "/" << max_runs << "... ";
            cout.flush();

            while(eval_count < max_evals) {
                for(int k=0; k<num_ants; ++k) {
                    if(eval_count >= max_evals) break;

                    construct_solution(ants[k]);
                    eval_count++;

                    if(ants[k].tour_length < run_best_dist) {
                        run_best_dist = ants[k].tour_length;
                    }
                    if(ants[k].tour_length < global_best_dist) {
                        global_best_dist = ants[k].tour_length;
                        global_best_path = ants[k].path;
                    }
                }
                update_pheromones(ants);
            }
            sum_best_dists += run_best_dist;
            cout << "Best: " << run_best_dist << endl;
        }
        
        // 計算 Mean Distance
        mean_best_dist = sum_best_dists / max_runs;

        cout << "------------------------------------------------" << endl;
        cout << "Average Best Distance: " << mean_best_dist << endl;
        cout << "Global Best Distance: " << global_best_dist << endl;
    }

    double get_best_dist() const { return global_best_dist; }
    double get_mean_dist() const { return mean_best_dist; } // 新增 Getter
    vector<int> get_best_path() const { return global_best_path; }
};

int main(int argc, char* argv[]) {
    string inputFileName;
    string outputFileName;
    
    if (argc > 2) {
        inputFileName = argv[1];
        outputFileName = argv[2];
    } else {
        cout << "Please input the input file name : ";
        cin >> inputFileName;
        cout << "Please input the output file name : ";
        cin >> outputFileName;
    }

    ifstream in(inputFileName);
    if (!in) {
        cerr << "Error: Cannot open input file." << endl;
        return 1;
    }

    vector<Point> cities;
    int pointNumber, x, y;
    while(in >> pointNumber >> x >> y) {
        cities.emplace_back(pointNumber, x, y);
    }
    in.close();

    if (cities.empty()) return 0;

    ACO aco(cities);
    aco.run();

    ofstream out(outputFileName);
    if (!out) return 1;

    // --- 修改：依照圖片格式輸出 Mean Distance 與 Best Distance ---
    out << "mean distance: " << fixed << setprecision(3) << aco.get_mean_dist() << endl;
    out << "distance: " << fixed << setprecision(3) << aco.get_best_dist() << endl;
    
    vector<int> path_indices = aco.get_best_path();
    for(int idx : path_indices) {
        out << cities[idx].PointNumber << endl;
    }
    
    out.close();
    cout << "Result saved to " << outputFileName << endl;

    return 0;
}