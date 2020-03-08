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
    double sum = 0;
    unsigned int seed = omp_get_wtime();
    int sign;   
    #pragma omp parallel for shared(A) firstprivate(seed)
    for (i = 0; i < leng; i++) {
        seed = omp_get_wtime();
        sign = 1;
        if ((rand_r(&seed) / (float) RAND_MAX) < 0.5) sign = -1;
        if ((wpart / (float) RAND_MAX) < 0.5) sign = -1;
        A[i].real(sign * (rand_r(&seed) / (float) RAND_MAX));
        A[i].imag(sign * (rand_r(&seed) / (float) RAND_MAX));
        sum += abs(A[i] * A[i]);
    }
    sum = sqrt(sum);
    #pragma omp parallel for 
    for (i = 0; i < leng; i++) {
        A[i] = A[i] / sum;
    }
    return A;
}

complexd* qubit_transform(complexd* A, int n, complexd* H, int k) {
    unsigned long long i, leng = 1 << n, ik = 1 << (n - k);
    complexd *B = new complexd[leng];
    #pragma omp parallel for shared(B)
        for (i = 0; i < leng; i++) {
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
    if (n < k) return 1;
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