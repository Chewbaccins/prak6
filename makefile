all: main.cpp
	g++ main.cpp -Wall -Werror -fopenmp -o run

q20_1: all
	echo "" >> results.txt
	echo "1 " >> results.txt
	for k in 1 10 20; do \
		for i in 1; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 20 $$k; \
			done \
		done \
	done

q20_2: all
	echo "" >> results.txt
	echo "2 " >> results.txt
	for k in 1 10 20; do \
		for i in 2; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 20 $$k; \
			done \
		done \
	done

q20_4: all
	echo "" >> results.txt
	echo "4 " >> results.txt
	for k in 1 10 20; do \
		for i in 4; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 20 $$k; \
			done \
		done \
	done

q20_8: all
	echo "" >> results.txt
	echo "8 " >> results.txt
	for k in 1 10 20; do \
		for i in 8; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 20 $$k; \
			done \
		done \
	done

q24_1: all
	echo "" >> results.txt
	echo "1 " >> results.txt
	for k in 1 10 24; do \
		for i in 1; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 24 $$k; \
			done \
		done \
	done

q24_2: all
	echo "" >> results.txt
	echo "2 " >> results.txt
	for k in 1 10 24; do \
		for i in 2; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 24 $$k; \
			done \
		done \
	done

q24_4: all
	echo "" >> results.txt
	echo "4 " >> results.txt
	for k in 1 10 24; do \
		for i in 4; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 24 $$k; \
			done \
		done \
	done

q24_8: all
	echo "" >> results.txt
	echo "8 " >> results.txt
	for k in 1 10 24; do \
		for i in 8; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 24 $$k; \
			done \
		done \
	done

q28_1: all
	echo "" >> results.txt
	echo "1 " >> results.txt
	for k in 1 10 28; do \
		for i in 1; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 28 $$k; \
			done \
		done \
	done

q28_2: all
	echo "" >> results.txt
	echo "2 " >> results.txt
	for k in 1 10 28; do \
		for i in 2; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 28 $$k; \
			done \
		done \
	done

q28_4: all
	echo "" >> results.txt
	echo "4 " >> results.txt
	for k in 1 10 28; do \
		for i in 4; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 28 $$k; \
			done \
		done \
	done

q28_8: all
	echo "" >> results.txt
	echo "8 " >> results.txt
	for k in 1 10 28; do \
		for i in 8; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 28 $$k; \
			done \
		done \
	done

q30_1: all
	echo "" >> results.txt
	echo "1 " >> results.txt
	for k in 1 10 30; do \
		for i in 1; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 30 $$k; \
			done \
		done \
	done

q30_2: all
	echo "" >> results.txt
	echo "2 " >> results.txt
	for k in 1 10 30; do \
		for i in 2; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 30 $$k; \
			done \
		done \
	done

q30_4: all
	echo "" >> results.txt
	echo "4 " >> results.txt
	for k in 1 10 30; do \
		for i in 4; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 30 $$k; \
			done \
		done \
	done

q30_8: all
	echo "" >> results.txt
	echo "8 " >> results.txt
	for k in 1 10 30; do \
		for i in 8; do \
			for n in 0 1 2 3 4 5 6 7 8 9; do \
				bsub -n 1 -W 00:25 -J kvant OMP_NUM_THREADS=$$i ./run 30 $$k; \
			done \
		done \
	done