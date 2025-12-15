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
// 用於產生 0.0 到 1.0 之間的隨機小數，用在輪盤賭選擇
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0.0, 1.0);

// --- Point 類別 (城市資料結構) ---
class Point {
public:
    int PointNumber; // 城市編號
    int X;           // X座標
    int Y;           // Y座標
    
    Point(int pointNumber = 0, int x = 0, int y = 0) {
        PointNumber = pointNumber;
        X = x;
        Y = y;
    }
};

// --- Ant 結構 (螞蟻個體) ---
struct Ant {
    vector<int> path;       // 螞蟻走過的路徑順序 (存城市 index)
    vector<bool> visited;   // 紀錄表：visited[i] = true 代表第 i 個城市去過了
    double tour_length;     // 這隻螞蟻走完一圈的總距離

    Ant(int n) {
        visited.resize(n, false);
        tour_length = 0.0;
        path.reserve(n + 1);
    }

    // 每回合開始前，重置螞蟻狀態
    void clear(int n) {
        path.clear();
        fill(visited.begin(), visited.end(), false);
        tour_length = 0.0;
    }
};

class ACO {
private:
    // ACO 參數
    int num_ants;       // 螞蟻數量
    double alpha;       // 費洛蒙權重 (History)
    double beta;        // 距離權重 (Heuristic)
    double rho;         // 費洛蒙揮發率
    double Q;           // 費洛蒙釋放常數
    int max_runs;       // 實驗總執行次數 (為了算平均表現)
    
    vector<Point> cities; // 城市列表
    int n;                // 城市總數
    
    vector<vector<double>> dist_matrix;     // 距離矩陣 (存任兩點距離)
    vector<vector<double>> pheromone;       // 費洛蒙矩陣 (存任兩點路徑上的氣味濃度)
    vector<vector<double>> heuristic_cache; // 啟發值矩陣 (預先算好 1/distance，加速運算)

    // 最佳解紀錄
    vector<int> global_best_path; // 目前找到最好的路徑
    double global_best_dist;      // 目前找到最短的距離
    
    // 平均距離紀錄 (作業要求)
    double mean_best_dist; 

public:
    // 建構子：初始化參數與矩陣
    ACO(const vector<Point>& pts) : cities(pts), n(pts.size()) {
        num_ants = 30;      // 每一輪派 30 隻螞蟻
        alpha = 1.0;        // 費洛蒙影響力
        beta = 2.0;         // 距離影響力 (傾向走近的)
        rho = 0.1;          // 每次揮發 10%
        Q = 100.0;          // 費洛蒙常數
        max_runs = 30;      // 跑 30 次實驗取平均
        
        global_best_dist = 1e18; // 設一個超大值，確保第一次就能更新
        mean_best_dist = 0.0;

        // 1. 預先計算距離矩陣 (避免重複算 sqrt，浪費時間)
        dist_matrix.resize(n, vector<double>(n));
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                double dx = (double)cities[i].X - cities[j].X;
                double dy = (double)cities[i].Y - cities[j].Y;
                dist_matrix[i][j] = sqrt(dx*dx + dy*dy);
            }
        }

        // 2. 預先計算啟發值 (Heuristic Information)
        // 公式： eta = 1.0 / distance
        // 意義：距離越短，值越大，代表螞蟻越想選
        heuristic_cache.resize(n, vector<double>(n));
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                if(i != j) {
                    double dist = dist_matrix[i][j];
                    if (dist < 1e-10) dist = 1e-10; // 避免除以 0
                    double eta = 1.0 / dist;
                    // 預先算好 eta^beta，之後選路時直接查表就好
                    heuristic_cache[i][j] = pow(eta, beta);
                }
            }
        }
    }

    // 初始化費洛蒙：每一回合開始，地圖都要是乾淨的(或設為初始值)
    void init_pheromones() {
        pheromone.assign(n, vector<double>(n, 1.0)); // 初始濃度設為 1.0
    }

    // 讓一隻螞蟻從起點走到終點，建構完整路徑
    void construct_solution(Ant& ant) {
        ant.clear(n);
        int start_node = rand() % n; // 隨機選一個起點
        ant.path.push_back(start_node);
        ant.visited[start_node] = true;

        int current = start_node;
        // 走 n-1 步，訪問完剩下所有城市
        for(int step = 0; step < n - 1; ++step) {
            int next_node = select_next_node(ant, current);
            ant.path.push_back(next_node);
            ant.visited[next_node] = true;
            ant.tour_length += dist_matrix[current][next_node];
            current = next_node;
        }
        // 最後一步：回到起點 (TSP 要求形成迴路)
        ant.tour_length += dist_matrix[current][start_node];
    }

    // 核心邏輯：輪盤賭選擇下一個城市
    int select_next_node(const Ant& ant, int current) {
        double sum_prob = 0.0;
        
        // 1. 計算分母：所有「未訪問」城市的吸引力總和
        for(int j=0; j<n; ++j) {
            if(!ant.visited[j]) { 
                double tau = pheromone[current][j];        // 這條路上的費洛蒙
                double eta_pow_beta = heuristic_cache[current][j]; // 這條路的距離啟發值
                
                // 吸引力公式： P = (費洛蒙^alpha) * (距離分數^beta)
                double p = (alpha == 1.0) ? (tau * eta_pow_beta) : (pow(tau, alpha) * eta_pow_beta);
                sum_prob += p;
            }
        }

        // 2. 輪盤賭 (Roulette Wheel Selection)
        // 概念：機率越大的城市，佔的圓餅圖面積越大，越容易被丟中的球 (r) 選到
        double r = dis(gen) * sum_prob; 
        double cum_prob = 0.0;
        
        for(int j=0; j<n; ++j) {
            if(!ant.visited[j]) {
                double tau = pheromone[current][j];
                double eta_pow_beta = heuristic_cache[current][j];
                double p = (alpha == 1.0) ? (tau * eta_pow_beta) : (pow(tau, alpha) * eta_pow_beta);
            
                cum_prob += p; // 累加機率
                if(r <= cum_prob) return j; // 如果球落在這個區間，就選這個城市
            }
        }
        return -1; // 理論上不該執行到這
    }

    // 更新費洛蒙：揮發 + 堆積
    void update_pheromones(const vector<Ant>& ants) {
        // 1. 全局揮發 (Evaporation)
        // 讓舊的費洛蒙隨時間變淡，避免收斂到局部最佳解
        double evaporation = 1.0 - rho;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                pheromone[i][j] *= evaporation; 
            }
        }

        // 2. 費洛蒙堆積 (Deposit)
        // 所有螞蟻走完後，根據它們的表現，在路徑上留下味道
        for(const auto& ant : ants) {
            // 路徑越短 (tour_length 小)，留下的量 (delta) 越多
            double delta = Q / ant.tour_length;
            
            // 沿著螞蟻走過的路徑加強費洛蒙
            for(size_t i=0; i < ant.path.size() - 1; ++i) {
                int u = ant.path[i];
                int v = ant.path[i+1];
                pheromone[u][v] += delta; 
                pheromone[v][u] += delta; // 雙向都要加，因為是無向圖
            }
            // 處理最後回到起點的那一段
            int u = ant.path.back();
            int v = ant.path[0];
            pheromone[u][v] += delta;
            pheromone[v][u] += delta;
        }
    }

    // 主執行迴圈
    void run() {
        cout << "Dataset size: " << n << " cities." << endl;
        // 計算總評估次數 (作業可能限制運算量)
        cout << "Max Evaluations per run: " << 10000 * n << endl;

        double sum_best_dists = 0; // 用來算平均值
        int max_evals = 10000 * n;
        vector<Ant> ants(num_ants, Ant(n));

        // 執行多次實驗 (max_runs = 30)
        for(int run = 0; run < max_runs; ++run) {
            init_pheromones(); // 每次實驗重新初始化費洛蒙
            double run_best_dist = 1e18; // 該次實驗的最佳解
            int eval_count = 0;

            cout << "Running " << setw(2) << run + 1 << "/" << max_runs << "... ";
            cout.flush();

            // 在運算次數限制內，不斷迭代
            while(eval_count < max_evals) {
                // 1. 讓所有螞蟻走一遍
                for(int k=0; k<num_ants; ++k) {
                    if(eval_count >= max_evals) break;

                    construct_solution(ants[k]); // 建構路徑
                    eval_count++;

                    // 更新該次實驗的最佳解
                    if(ants[k].tour_length < run_best_dist) {
                        run_best_dist = ants[k].tour_length;
                    }
                    // 更新全域歷史最佳解
                    if(ants[k].tour_length < global_best_dist) {
                        global_best_dist = ants[k].tour_length;
                        global_best_path = ants[k].path;
                    }
                }
                // 2. 所有螞蟻走完後，更新環境費洛蒙
                update_pheromones(ants);
            }
            sum_best_dists += run_best_dist;
            cout << "Best: " << run_best_dist << endl;
        }
        
        // 計算 30 次實驗的平均最佳距離
        mean_best_dist = sum_best_dists / max_runs;

        cout << "------------------------------------------------" << endl;
        cout << "Average Best Distance: " << mean_best_dist << endl;
        cout << "Global Best Distance: " << global_best_dist << endl;
    }

    // Getters，用於輸出結果
    double get_best_dist() const { return global_best_dist; }
    double get_mean_dist() const { return mean_best_dist; } 
    vector<int> get_best_path() const { return global_best_path; }
};

int main(int argc, char* argv[]) {
    // 處理輸入輸出檔名
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

    // 讀取檔案
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

    // 執行 ACO 演算法
    ACO aco(cities);
    aco.run();

    // 輸出結果檔案
    ofstream out(outputFileName);
    if (!out) return 1;

    // 輸出格式：先平均距離，再最佳距離，最後是路徑
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