#include <iostream>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include "omp.h"
#include "time.h"
#include "sys/time.h"

#define eps 0.01
using namespace std;

typedef complex<double> complexd;

complexd *generate(unsigned long long seg_leng, int rank) {
    auto *gen_buf = new complexd[seg_leng];
    double sqr = 0, module;
    unsigned int seed = rand() + rank;
#pragma omp parallel for schedule(static) shared(gen_buf) reduction(+: sqr)
        for (unsigned long long i = 0; i < seg_leng; i++) {
            gen_buf[i].real((rand_r(&seed) / (float) RAND_MAX) - 0.5f);
            gen_buf[i].imag((rand_r(&seed) / (float) RAND_MAX) - 0.5f);
            sqr += abs(gen_buf[i] * gen_buf[i]);
        }
    MPI_Reduce(&sqr, &module, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        module = sqrt(module);
    }
    MPI_Bcast(&module, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#pragma omp parallel for schedule(static) shared(gen_buf, module)
        for (std::size_t i = 0; i < seg_leng; i++) {
            gen_buf[i] /= module;
        }
    return gen_buf;
}

complexd *read(char *f, unsigned long long seg_leng, int rank) {
    MPI_File file;
    if (MPI_File_open(MPI_COMM_WORLD, f, MPI_MODE_RDONLY, MPI_INFO_NULL, &file)) {
        if (!rank)
            printf("Error opening file %s\n", f);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    auto *proc_buf = new complexd[seg_leng];
    double read_buf[2];
    MPI_File_seek(file, 2 * seg_leng * rank * sizeof(double), MPI_SEEK_SET);
    for (std::size_t i = 0; i < seg_leng; ++i) {
        MPI_File_read(file, &read_buf, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
        proc_buf[i].real(read_buf[0]);
        proc_buf[i].imag(read_buf[1]);
    }
    MPI_File_close(&file);
    return proc_buf;
}


void write(char *output, complexd *out_buf, unsigned long long seg_leng, int rank) {
    MPI_File file;
    if (MPI_File_open(MPI_COMM_WORLD, output, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                      MPI_INFO_NULL, &file)) {
        if (!rank)
            printf("Error opening file %s\n", output);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    double proc_buf[2];
    MPI_File_seek(file, 2 * seg_leng * rank * sizeof(double), MPI_SEEK_SET);
    for (std::size_t i = 0; i < seg_leng; ++i) {
        proc_buf[0] = out_buf[i].real();
        proc_buf[1] = out_buf[i].imag();
        MPI_File_write(file, &proc_buf, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&file);
}

void
qubit_transform(complexd *buf_zone, complexd U[2][2], unsigned long long seg_leng, unsigned int k, 
                complexd *recv_zone, int rank) {
    unsigned first_index = rank * seg_leng;
    int rank_neighbor = first_index ^(1u << (k - 1));
    rank_neighbor /= seg_leng;
    if (rank != rank_neighbor) {
        MPI_Sendrecv(buf_zone, seg_leng, MPI_DOUBLE_COMPLEX, rank_neighbor, 0, recv_zone, seg_leng, MPI_DOUBLE_COMPLEX,
                    rank_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank > rank_neighbor) {
#pragma omp parallel for schedule(static) shared(recv_zone, buf_zone, U)
            for (unsigned long long i = 0; i < seg_leng; i++) {
                recv_zone[i] = U[1][0] * recv_zone[i] + U[1][1] * buf_zone[i];
            }
        } else {
#pragma omp parallel for schedule(static) shared(recv_zone, buf_zone, U)
            for (unsigned long long i = 0; i < seg_leng; i++) {
                recv_zone[i] = U[0][0] * buf_zone[i] + U[0][1] * recv_zone[i];
            }
        }
    } else {
        unsigned shift = (int) log2(seg_leng) - k;
        unsigned pow = 1u << (shift);
#pragma omp parallel for schedule(static) shared(recv_zone, buf_zone, U)
        for (std::size_t i = 0; i < seg_leng; i++) {
            unsigned i0 = i & ~pow;
            unsigned i1 = i | pow;
            unsigned iq = (i & pow) >> shift;
            recv_zone[i] = U[iq][0] * buf_zone[i0] + U[iq][1] * buf_zone[i1];
        }
    }
}

double normal_dis_gen()
{
    double S = 0.;
    for (int i = 0; i < 12; ++i) { S += (double) rand() / RAND_MAX; }
    return S - 6.0;
}

void add_noise(complexd U[2][2], complexd V[2][2], double tet) {
    V[0][0] = U[0][0] * cos(tet) - U[0][1] * sin(tet);
    V[0][1] = U[0][0] * sin(tet) + U[0][1] * cos(tet);
    V[1][0] = U[1][0] * cos(tet) - U[1][1] * sin(tet);
    V[1][1] = U[1][0] * sin(tet) + U[1][1] * cos(tet);
}

double one_f(complexd *ideal, complexd *noise, int rank, unsigned long long seg_size) {
    double sqr = 0;
    double module;
    for (unsigned long long i = 0; i < seg_size; i++) {
        sqr += abs(ideal[i] * conj(noise[i])) * abs(ideal[i] * conj(noise[i]));
    }
    MPI_Reduce(&sqr, &module, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        return module;
    } else {
        return 0;
    }
}

int main(int argc, char **argv) {
    //mode 0 - read file, 1 - generate data
    //example: mpirun ./a.out n UNSIGNED INPUT_FILE OUTPUT_FILE
    //or
    //example: mpirun ./a.out n UNSIGNED OUTPUT_FILE
    bool file_input = false;
    char *input, *output;
    unsigned n;
    n = atoi(argv[2]);
    if (file_input == true) {
        input = argv[3];
        output = argv[4];
    } else {
        output = argv[3]; 
    }

    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    unsigned long long whole_leng = 1LLU << n;
    unsigned long long seg_leng = whole_leng / size;
    complexd *input_buf;
    auto *norm_buf = new complexd[seg_leng];
    auto *noise_buf = new complexd[seg_leng];
    if (!file_input) {
        input_buf = generate(seg_leng, rank);
    } else {
        input_buf = read(input, seg_leng, rank);
    }
    
    complexd U[2][2];
    U[0][0] = 1 / sqrt(2);
    U[0][1] = 1 / sqrt(2);
    U[1][0] = 1 / sqrt(2);
    U[1][1] = -1 / sqrt(2);

    //main
    double begin = MPI_Wtime();
    for (unsigned qubit = 1; qubit < n + 1; qubit++) {
        qubit_transform(input_buf, U, seg_leng, qubit, norm_buf, rank);
    }
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Main time " << end - begin << " seconds\n";
    }

    //noise
    double noise_time = 0;
    for (unsigned qubit = 1; qubit < n + 1; qubit++) {
        complexd U_changed[2][2];
        double tet = 0;
        if (rank == 0) {
            tet = normal_dis_gen() * eps;
        }
        MPI_Bcast(&tet, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        add_noise(U, U_changed, tet);
        begin = MPI_Wtime();
        qubit_transform(input_buf, U, seg_leng, qubit, noise_buf, rank);
        end = MPI_Wtime();
        noise_time += end - begin;
    }
    if (rank == 0) {
        std::cout << "Noise time " << noise_time << " seconds\n";
    }
    double distance = one_f(norm_buf, noise_buf, rank, seg_leng);
    if (!rank) {
        cout << "distance " << distance << endl;
    }

    MPI_Finalize();
    delete[] input_buf;
}