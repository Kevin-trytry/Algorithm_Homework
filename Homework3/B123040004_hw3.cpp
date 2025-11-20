#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<climits>
#include<algorithm>
#include<iomanip> // 為了 setprecision

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
        TSP() : MinDistance(1e18) {} // 使用更大的初始值

        void SolveDP(const vector<Point> &cities) {
            if (cities.size() <= 1) {
                MinDistance = (cities.size() == 1) ? 0.0 : 0.0;
                MinDistanceResult = cities;
                if (cities.size() == 1) MinDistanceResult.push_back(cities[0]); // 閉合迴路
                return;
            }

            int n = cities.size();
            // dist[i][j] 儲存城市 i 到城市 j 的距離
            vector<vector<double>> dist(n, vector<double>(n));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    dist[i][j] = sqrt(pow((cities[i].X - cities[j].X), 2) + pow((cities[i].Y - cities[j].Y), 2));
                }
            }

            // DP table: dp[mask][i]
            vector<vector<double>> dp(1 << n, vector<double>(n, MinDistance));
            // Parent table: parent[mask][i] 儲存到達狀態 (mask, i) 的前一個城市
            vector<vector<int>> parent(1 << n, vector<int>(n, -1));

            // 1. 初始化：從城市 0 開始
            dp[1][0] = 0.0;

            // 2. 迭代計算 DP 表
            for (int mask = 1; mask < (1 << n); ++mask) {
                for (int j = 0; j < n; ++j) {
                    // 如果城市 j 在 mask 中
                    if (mask & (1 << j)) {
                        int prev_mask = mask ^ (1 << j); // 從 mask 中移除 j
                        
                        // 檢查 prev_mask 是否為空 (只剩 j) 或 prev_mask 包含起始城市 0
                        if (prev_mask == 0) continue; 

                        // 遍歷所有可能的前一個城市 i
                        for (int i = 0; i < n; ++i) {
                            // 城市 i 在 prev_mask 中
                            if (i != j && (prev_mask & (1 << i))) {
                                // 狀態轉移
                                double new_dist = dp[prev_mask][i] + dist[i][j];
                                if (new_dist < dp[mask][j]) {
                                    dp[mask][j] = new_dist;
                                    parent[mask][j] = i; // 記錄最佳前驅城市
                                }
                            }
                        }
                    }
                }
            }

            // 3. 找出最短距離與最佳終點
            int final_mask = (1 << n) - 1;
            int best_last_city = -1;
            MinDistance = 1e18;

            for (int j = 1; j < n; ++j) { // 終點城市 j 不能是起始城市 0
                double final_dist = dp[final_mask][j] + dist[j][0];
                if (final_dist < MinDistance) {
                    MinDistance = final_dist;
                    best_last_city = j;
                }
            }

            // 4. 路徑回溯
            MinDistanceResult.clear();

            if (best_last_city != -1) {
                // 路徑重建：從終點 best_last_city 逆向回溯到城市 0
                vector<int> path_indices;
                int current_mask = final_mask;
                int current_city = best_last_city;

                while (current_city != -1) {
                    path_indices.push_back(current_city);
                    int prev_city = parent[current_mask][current_city];
                    
                    if (prev_city == -1 && current_city == 0 && current_mask == 1) {
                        // 達到起始狀態 (mask=1, city=0)
                        break;
                    } else if (prev_city == -1) {
                        // 不可達的狀態，可能計算有誤，但理論上不該發生
                        break; 
                    }
                    
                    current_mask = current_mask ^ (1 << current_city);
                    current_city = prev_city;
                }

                // path_indices 現在是 (Pn-1, Pn-2, ..., P1, 0)
                // 閉合迴路是： 0 -> P1 -> ... -> Pn-1 -> 0
                
                // 1. 起始城市 (0)
                for(int i = path_indices.size() - 1; i >= 0; --i) {
                    MinDistanceResult.push_back(cities[path_indices[i]]);
                }
                
                // 2. 回到起始城市 (0) 完成閉合
                MinDistanceResult.push_back(cities[0]);
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

// 保持 main 函數不變 (增加檔案檢查)
int main( ) {
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

    // 檢查輸入檔是否開啟正確
    if (in.fail()) {
        cout << "Error: Input file \"" << inputFileName << "\" cannot open or does not exist." << endl;
        return 1;
    }

    int x, y, pointNumber;
    while(in >> pointNumber >> x >> y) {
        // 重要：cities 向量的索引是 0 到 N-1，因此 PointNumber 應該從 0 開始連續編號，
        // 否則 DP 邏輯會錯亂。如果您的輸入 PointNumber 不是從 0 開始，這裡需要轉換。
        // 但由於原始 DP 邏輯是依賴城市在 cities 向量中的索引，我們維持不變。
        Point readin(cities.size(), x, y); // 這裡改為使用 cities.size() 作為 DP 內部索引
        readin.PointNumber = pointNumber;  // 保留原始的 PointNumber
        cities.push_back(readin);
    }
    in.close();

    // 檢查是否成功讀取到城市資料
    if (cities.empty()) {
        cout << "Error: No cities were read from the input file. Check file format." << endl;
        return 1;
    }

    // implement tsp using DP
    TSP tsp;
    tsp.SolveDP(cities); // 呼叫 DP 求解函數

    double minDistance = tsp.GetMinDistance();
    vector<Point> minPermResult = tsp.GetMinDistanceResult();

    // output to the ans.txt
    out.open(outputFileName);

    // 檢查輸出檔是否開啟正確
    if (out.fail()) {
        cout << "Error: Output file \"" << outputFileName << "\" cannot be created or opened." << endl;
        return 1;
    }

    // 格式化輸出：總距離 (保留小數點)
    out << fixed << setprecision(10) << minDistance << endl;

    // 輸出路徑
    for (const auto &var : minPermResult)
        out << var.PointNumber << endl;
    out.close();
    
    cout << "Successfully calculated TSP. Result written to " << outputFileName << endl;

    return 0;
}