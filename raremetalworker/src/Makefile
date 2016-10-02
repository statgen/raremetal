EXE=raremetalworkerTest
#OMP_EXE= raremetalworker-omp
TOOLBASE = FastFit TransformResiduals KinshipEmp OutputKin PreMeta AutoFit AutoFit2 DataQC WritePDF GroupFromAnnotation logistic
SRCONLY = Main.cpp 
HDRONLY = WriteLog.h snpHWE.h

########################
# Handle openmp
ifeq ($(MAKECMDGOALS),openmp)
  OBJDIR_OPT = $(OBJDIR)/omp
#  EXE = $(OMP_EXE)
endif
########################

USER_INCLUDES = -I ../../otherLib/eigen-3.2.0 -I$(LIB_RAREMETAL_PATH)include
LIB_RAREMETAL_PATH = ../../libRareMetal/
RMATH_LIB =  $(LIB_RAREMETAL_PATH)/libRmath.a

# USER_WARNINGS ?= -Werror -Wno-unused-variable $(shell if [ X$(CCVERSION) \> X4.2.0 ] ; then echo " -Wno-unused-result" ; fi)
USER_WARNINGS =

COMPILE_ANY_CHANGE = 

# Set additional values prior to including the base Makefiles.
# Set dependencies
opt: $(LIB_RAREMETAL_PATH)libRareMetal.a 
debug: $(LIB_RAREMETAL_PATH)libRareMetal_debug.a
profile: $(LIB_RAREMETAL_PATH)libRareMetal_profile.a

# Define USER_LIBS as it should be used in each build type
opt: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal.a $(RMATH_LIB)
debug: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal_debug.a $(RMATH_LIB)
profile: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal_profile.a $(RMATH_LIB)

PARENT_MAKE = Makefile.src
include ../../Makefile.inc

########################
# Add any additional targets here.

$(LIB_RAREMETAL_PATH)libRareMetal.a:
	echo here
	$(MAKE) -C $(@D) --no-print-directory opt
$(LIB_RAREMETAL_PATH)libRareMetal_debug.a:
	echo here
	$(MAKE) -C $(@D) --no-print-directory debug
$(LIB_RAREMETAL_PATH)libRareMetal_profile.a:
	$(MAKE) -C $(@D) --no-print-directory profile

USER_COMPILE_VARS = -DVERSION="\"${VERSION}\"" -g -O0 -fstack-protector-all

#binary: USER_COMPILE_VARS += -static
#binary: USER_COMPILE_VARS += -lstdc++
binary: opt

openmp: USER_COMPILE_VARS += -fopenmp
openmp: opt debug profile

USER_REMOVES += -rm -rf $(OBJDIR)/omp/*.o
