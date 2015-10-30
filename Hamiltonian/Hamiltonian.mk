LSPATH := ./Hamiltonian

LLIBS := VCA_Fortran
LIBS += $(LLIBS)
$(LLIBS):
	$(FCC) -c $(LSPATH)/$@.f90 -m $@
	mv *.so $(LSPATH)/
