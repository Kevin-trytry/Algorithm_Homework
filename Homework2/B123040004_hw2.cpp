#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<climits>
#include<algorithm>
using namespace std;

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

// Greedy Algorithm solving for TSP
class TSP {
    public:
        TSP() {
            MinDistance = INT_MAX;
        }
        void Greedy(vector<Point> &cities, vector<Point> &visited, double distance, Point StartPoint) {
            // Link the last point to the start point
            if (visited.size() == cities.size()) {
                distance += Distance(visited.back(), StartPoint);
                if (distance < MinDistance) {
                    MinDistance = distance;
                    MinDistanceResult = visited;
                }
                return;
            }
            // Greedy to find the nearest point
            double two_point_mindistance = INT_MAX;
            Point nearest_point;
            for (auto &var : cities) {
                if (Find(var, visited)) continue;
                double newDistance = Distance(visited.back(), var);
                if (newDistance < two_point_mindistance) {
                    two_point_mindistance = newDistance;
                    nearest_point = var;   
                }
            }
            visited.push_back(nearest_point);
            Greedy(cities, visited, distance + two_point_mindistance, StartPoint);
        }

        /// @brief 
        /// @param point 
        /// @param visited 
        /// @return if the point is in the visited vector
        bool Find(Point point, vector<Point> &visited) {
            for (auto &var : visited) {
                if (var.PointNumber == point.PointNumber) return true;
            }
            return false;
        }

        double Distance(Point &a, Point &b) {
            return sqrt(pow((a.X - b.X), 2) + pow((a.Y - b.Y), 2));
        }

        double GetMinDistance() {
            return MinDistance;
        }   
        vector<Point> GetMinDistanceResult() {
            return MinDistanceResult;
        }
    private:
        double MinDistance;
        vector<Point> MinDistanceResult;
        vector<Point> visited;
};

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

    // check if the file open correctly
    if (in.fail()) cout << "File cannot open." << endl;

    int x, y, pointNumber;
    while(in >> pointNumber >> x >> y) {
        Point readin(pointNumber, x, y);
        cities.push_back(readin);
    }
    in.close();

    // implement tsp
    TSP tsp;
    // try each point as the start point
    for (auto &var : cities) {
        vector<Point> visited;
        visited.push_back(var);
        tsp.Greedy(cities, visited, 0, var);
    }
    double minDistance = tsp.GetMinDistance();
    vector<Point> minPermResult = tsp.GetMinDistanceResult();
    
    //output to the ans.txt
    out.open(outputFileName);
    out << minDistance << endl;
    for (auto &var : minPermResult) 
        out << var.PointNumber << endl;
    out.close();
    
    return 0;
}
