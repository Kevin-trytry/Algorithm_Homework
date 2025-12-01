import matplotlib.pyplot as plt
import numpy as np
import sys

def read_ans(filename):
    ans = []
    # Reading answers from the file
    try :
        with open(filename, 'r') as f:
            for line in f.readlines()[1:]:
                ans.extend(list(map(int, line.split())))
    except:
        print("Error: Cannot open the answer file or file does not exist.")
        sys.exit(1)
    return ans

def read_points(filename):
    x_coords = []
    y_coords = []
    # Reading points from the file
    try:
        with open(filename, 'r') as f:
            for line in f.readlines():
                pointNumber, x, y = map(float, line.strip().split())
                pointNumber = int(pointNumber)
                x_coords.append(x)
                y_coords.append(y)
    except:
        print("Error: Cannot open the point file or file does not exist.")
        sys.exit(1)
    return x_coords, y_coords

def plot_png(point_file, ans_file, fig_file):
    # Reading points and answers
    x_coords, y_coords = read_points(point_file)
    ans = read_ans(ans_file)
    

    # Draw TSP path
    path_x = [x_coords[i - 1] for i in ans] + [x_coords[ans[0] - 1]]
    path_y = [y_coords[i - 1] for i in ans] + [y_coords[ans[0] - 1]]
    plt.plot(path_x, path_y, '-', color='blue', linewidth=1.5, label='TSP Path')
    

    # Plotting the points
    plt.xlim(0, 100) 
    plt.ylim(0, 100)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('TSP Solution Visualization')
    plt.grid(True)
    plt.savefig(fig_file, dpi=300)
    plt.show()

if (__name__ == "__main__"):
    point_file = input("Please input the point file name(e.g., point.txt): ")
    ans_file = input("Please input the answer file name(e.g., ans.txt): ")
    fig_file = input("Please input the figure file name you want to output(e.g., fig.png): ")
    plot_png(point_file, ans_file, fig_file)
