all:
	mpic++ main.cpp -fopenmp
	
test:
	for n in 20 24 28 30; do \
		for k in 1 3 $$n; do \
			for i in 1 2 4 8; do \
				bsub -W 15 -q normal OMP_NUM_THREADS=8  mpiexec -n $$i ./a.out n $$n k $$k out.txt \
			done \
		done \
	done