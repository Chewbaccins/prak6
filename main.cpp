#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <omp.h>
using namespace std;

typedef complex<double> complexd;

complexd* createq(int n) {
    unsigned long long leng = 1 << n, i;
    complexd *A = new complexd[leng];
    unsigned int seed = omp_get_wtime(), wpart;
    int sign;   
    #pragma omp parallel for shared(A, leng) private(i) firstprivate(seed)
    for (i = 0; i < leng; i++) {
        seed = omp_get_wtime();
        sign = 1;
        wpart = rand_r(&seed);
        if ((wpart / (float) RAND_MAX) < 0.5) {
            sign = -1;
        }
        A[i].real(sign * (wpart + rand_r(&seed) / (float) RAND_MAX));
        
        sign = 1;
        wpart = rand_r(&seed);
        if ((wpart / (float) RAND_MAX) < 0.5) {
            sign = -1;
        }
        A[i].imag(sign * (wpart + rand_r(&seed) / (float) RAND_MAX));
    }
    return A;
}

complexd* qubit_transform(complexd* A, int n, complexd* H, int k) {
    unsigned long long  leng = 1LLU << n, ik = 1LLU<< (n - k);
    complexd *B = new complexd[leng];
    #pragma omp parallel for shared(A, B, H, leng, ik)
        for (unsigned long long i = 0; i < leng; i++) {
            int num = i & ik;
            switch (num) {
            case 0:
                B[i] = H[0]*A[i] + H[1]*A[i | ik];
                break;
            default:
                B[i] = H[2]*A[i & ~ik] + H[3]*A[i];
                break;
            }
        }
    return B;
}

int main(int argc, char **argv) {

    if (argc < 3) {
        cout << "n and k required" << endl;
        return 0;
    }
    int n, k;
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    complexd *A = createq(n);
    complexd H[] = {1/sqrt(2), 1/sqrt(2), 1/sqrt(2), -1/sqrt(2)};

    double time = omp_get_wtime(); 
    complexd *B = qubit_transform(A, n, H, k);
    time = omp_get_wtime() - time;
    
    ofstream results;
    results.open ("results.txt", ios::out | ios::app);
    results << n << " " << k << " " << time << "\n";
    results.close();
    
    delete[] A;
    delete[] B;
    return 0;
}