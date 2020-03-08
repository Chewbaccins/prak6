all: main.cpp
	g++ main.cpp -Wall -Werror -fopenmp -O3 -o run

start: all
	for n in 20 24 28 30; do \
		for k in 1 10 $$n; do \
			for i in 1 2 4 8; do \
				for iter in 0 1 2 3 4; do \
					bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run $$n $$k; \
				done \
			done \
		done \
	done