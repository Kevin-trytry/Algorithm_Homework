#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<climits>
#include<algorithm>
#include<iomanip>

using namespace std;

// 保持 Point 類別不變
class Point {
public:
    Point(int pointNumber = 0, int x = 0, int y = 0) {
        PointNumber = pointNumber;
        X = x;
        Y = y;
    }
    int PointNumber;
    int X;
    int Y;
};

// --- 使用 Dynamic Programming (DP) 解決 TSP 問題 ---
class TSP {
public:
    TSP() : MinDistance(1e18) {}

    void SolveDP(const vector<Point> &cities) {
        int n = cities.size();
        if (n == 0) return;

        // 特殊情況處理：只有一個城市
        if (n == 1) {
            MinDistance = 0;
            MinDistanceResult = cities;
            return;
        }

        // 1. 預計算距離矩陣
        vector<vector<double>> dist(n, vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                dist[i][j] = sqrt(pow((double)cities[i].X - cities[j].X, 2) + 
                                  pow((double)cities[i].Y - cities[j].Y, 2));
            }
        }

        // 2. 初始化 DP 表
        // dp[mask][i] 表示經過 mask 集合中的城市，且最後停在城市 i 的最短距離
        // 1 << n 代表 2^n 種狀態
        vector<vector<double>> dp(1 << n, vector<double>(n, 1e18));
        
        // parent[mask][i] 用來回溯路徑，記錄到達該狀態的上一個城市
        vector<vector<int>> parent(1 << n, vector<int>(n, -1));

        // Base Case: 假設從城市 0 出發 (mask = 1, 二進位 00...01)
        dp[1][0] = 0;

        // 3. DP 迭代
        // mask 從 1 (00...01) 到 (1<<n)-1 (11...11)
        for (int mask = 1; mask < (1 << n); ++mask) {
            // 優化：如果 mask 不包含起始點 0，則此路徑無效 (因為我們固定從 0 出發)
            if (!(mask & 1)) continue;

            // i 是當前路徑的最後一個城市
            for (int i = 0; i < n; ++i) {
                // 如果 mask 包含城市 i
                if ((mask & (1 << i))) {
                    
                    // prev_mask 是還沒訪問 i 之前的狀態
                    int prev_mask = mask ^ (1 << i);
                    if (prev_mask == 0) continue; // 只有起點的情況已在 Base Case 處理

                    // j 是上一個城市
                    for (int j = 0; j < n; ++j) {
                        // 如果 prev_mask 包含 j
                        if ((prev_mask & (1 << j))) {
                            // 狀態轉移方程：嘗試從 j 走到 i
                            if (dp[prev_mask][j] + dist[j][i] < dp[mask][i]) {
                                dp[mask][i] = dp[prev_mask][j] + dist[j][i];
                                parent[mask][i] = j;
                            }
                        }
                    }
                }
            }
        }

        // 4. 計算回到起點的最短距離 (形成閉環)
        // full_mask 代表所有城市都已訪問 (11...11)
        int full_mask = (1 << n) - 1;
        MinDistance = 1e18;
        int last_city = -1;

        // 嘗試從任意最後一個城市 k 回到起點 0
        for (int k = 1; k < n; ++k) {
            double current_dist = dp[full_mask][k] + dist[k][0];
            if (current_dist < MinDistance) {
                MinDistance = current_dist;
                last_city = k;
            }
        }

        // 5. 路徑回溯
        MinDistanceResult.clear();
        if (last_city != -1) {
            vector<int> path_indices;
            int curr_mask = full_mask;
            int curr_city = last_city;

            while (curr_city != -1) {
                path_indices.push_back(curr_city);
                int temp_city = curr_city;
                curr_city = parent[curr_mask][curr_city];
                curr_mask = curr_mask ^ (1 << temp_city);
            }
            
            // 因為是從終點回溯，所以需要反轉
            reverse(path_indices.begin(), path_indices.end());

            // 將索引轉換回 Point 物件
            for (int idx : path_indices) {
                MinDistanceResult.push_back(cities[idx]);
            }
        }
    }

    double GetMinDistance() const {
        return MinDistance;
    }

    vector<Point> GetMinDistanceResult() const {
        return MinDistanceResult;
    }

private:
    double MinDistance;
    vector<Point> MinDistanceResult;
};

// --- Main 函數 (保留原始讀檔邏輯，僅微調輸出格式) ---
int main() {
    string inputFileName;
    string outputFileName;
    vector<Point> cities;

    cout << "Please input the input file name : ";
    cin >> inputFileName;
    cout << "Please input the output file name : ";
    cin >> outputFileName;

    // read in file
    ifstream in;
    ofstream out;
    in.open(inputFileName);

    if (in.fail()) {
        cout << "Error: Input file cannot open." << endl;
        return 1;
    }

    int x, y, pointNumber;
    while(in >> pointNumber >> x >> y) {
        Point readin(cities.size(), x, y);
        readin.PointNumber = pointNumber; 
        cities.push_back(readin);
    }
    in.close();

    // implement tsp using DP
    TSP tsp;
    tsp.SolveDP(cities);

    double minDistance = tsp.GetMinDistance();
    vector<Point> minPermResult = tsp.GetMinDistanceResult();

    // output to the ans.txt
    out.open(outputFileName);
    out << "distance: " << fixed << setprecision(3) << minDistance << endl;

    for (const auto &var : minPermResult)
        out << var.PointNumber << endl;
    
    out.close();

    cout << "Finished. Result saved to " << outputFileName << endl;

    return 0;
}