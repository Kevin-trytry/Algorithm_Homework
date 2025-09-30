#include <iostream>
using namespace std;

void Permutations (char *a, const int k, const int m)
//Generate all the permutations of a[k], ..., a[m]
{
    if (k == m) {  //Output permutation
     for (int i = 0; i <= m; i++) cout << a[i] << " ";
     cout << endl;
    }
    else { //a[k], ..., a[m] has more than one permutation
		for (int i = k; i <= m; i++) 
		{
		   swap(a[k], a[i]); // exchange
 		   Permutations(a, k+1, m);
		   swap(a[k], a[i]);
 		}
    } // end of else
} 

int main( ) {
    char test[4] = {'a', 'b', 'c', 'd'};
    Permutations(test, 0, 3);
}

