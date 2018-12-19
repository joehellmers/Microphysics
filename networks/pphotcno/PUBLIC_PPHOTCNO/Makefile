F90COMP = gfortran
SWITCHES = -ffree-line-length-none -O3 -fopenmp -Wall
OUTDIR = ~/bin
EXE = pphotcno

pphotcno: burn_common.dek const.dek helm_table_storage.dek implno.dek network.dek \
	public_pphotcno.f90 sparse_matrix.dek timers.dek vector_eos.dek
	$(F90COMP) public_pphotcno.f90 -o $(OUTDIR)/$(EXE) $(SWITCHES)

clean:
	rm -rf foo_000*.dat foo_summary.dat
	rm -rf avl_data.mod avl_stack.mod
	rm -rf $(OUTDIR)/$(EXE)
	rm -rf test_*
	rm -rf *.dat
