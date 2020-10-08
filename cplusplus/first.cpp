#include <iostream>

using namespace std;

// log2(3) = 1, log2(4) = 2, 
int log2(int i)
{
    int r = 0;
    while (i >>= 1) r++;
    return r;
}


int main(){
    int a = log2(10);
    cout << a << endl;    // cout needs to using namespace std
}