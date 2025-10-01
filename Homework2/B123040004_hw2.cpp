#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<climits>
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
            MinDistance = 1000000;
        }
        
        void Greedy(vector<Point> &cities, vector<Point> &visited, double distance, Point StartPoint) {
        
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