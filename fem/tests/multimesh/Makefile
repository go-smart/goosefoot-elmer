# A simple heat equation test
#
ELMER_GRID=/home/pmariappan/NUMABin/bin/ElmerGrid -metis 16
ELMER_SOLVER=/home/pmariappan/NUMABin/bin/ElmerSolver_mpi
run:
	$(ELMER_GRID) 1 2 global
	$(ELMER_GRID) 1 2 local
	$(ELMER_GRID) 1 2 reference
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out 
	/bin/rm -r global local reference
