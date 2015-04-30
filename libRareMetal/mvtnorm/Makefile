TOOLBASE = 
SRCONLY = randomF77.c mvt.f


PARENT_MAKE = Makefile.extlib
include ../Makefile.inc

$(OBJDIR_OPT)/%.o: %.f
	gfortran -c $*.f -o $@
$(OBJDIR_DEBUG)/%.o: %.f
	gfortran -c $*.f -o $@
$(OBJDIR_PROFILE)/%.o: %.f
	gfortran -c $*.f -o $@
