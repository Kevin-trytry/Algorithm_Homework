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

class TSP {
    public:
        TSP() {
            MinDistance = 1000000;
        }
        void Premutation(vector<Point> &cities, int start, int end) {
            if (start == end) 
                Distance(cities);
            else {
                for (int i = start; i <= end; i++) {
                    swap(cities[start], cities[i]);
                    Premutation(cities, start + 1, end);
                    swap(cities[start], cities[i]);
                }
            }
        }

        void Distance(vector<Point> &cities) {
            double totalDistance = 0;

            for (int i = 0; i < cities.size() - 1; i++) {
                totalDistance += sqrt(pow((cities[i].X - cities[i + 1].X), 2) + pow((cities[i].Y - cities[i + 1].Y), 2));
            }
            
            totalDistance += sqrt(pow((cities[0].X - cities[cities.size() - 1].X), 2) + pow((cities[0].Y - cities[cities.size() - 1].Y), 2));
            if (totalDistance < MinDistance) {
                MinDistance = totalDistance;
                MinDistanceResult = cities;
            }
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
    tsp.Premutation(cities, 0, cities.size() - 1);
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