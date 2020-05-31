#include <iostream>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include "omp.h"
#include "time.h"
#include "sys/time.h"

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

int main(int argc, char **argv) {
    //mode 0 - read file, 1 - generate data
    //example: mpirun a.out n UNSIGNED k UNSIGNED INPUT_FILE OUTPUT_FILE
    //example: mpirun a.out n UNSIGNED k UNSIGNED OUTPUT_FILE
    bool file_input = false;
    char *input, *output;
    unsigned k, n;
    n = atoi(argv[2]);
    k = atoi(argv[4]);
    if (file_input == true) {
        input = argv[5];
        output = argv[6];
    } else {
        output = argv[5]; 
    }

    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    unsigned long long whole_leng = 1LLU << n;
    unsigned long long seg_leng = whole_leng / size;
    complexd *input_buf;
    auto *recv_buf = new complexd[seg_leng];
    if (!file_input) {
        input_buf = generate(seg_leng, rank);
    } else {
        input_buf = read(input, seg_leng, rank);
    }
    write(output, input_buf, n, rank);
    
    complexd U[2][2];
    U[0][0] = 1 / sqrt(2);
    U[0][1] = 1 / sqrt(2);
    U[1][0] = 1 / sqrt(2);
    U[1][1] = -1 / sqrt(2);

    double begin = MPI_Wtime();
    qubit_transform(input_buf, U, seg_leng, k, recv_buf, rank);
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << "The process took " << end - begin << " seconds to run." << std::endl;
    }
    write(output, recv_buf, n, rank);
    MPI_Finalize();
    delete[] input_buf;
}