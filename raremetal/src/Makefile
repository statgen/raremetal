EXE=raremetal
TOOLBASE = Meta GroupFromAnnotation WritePDF SummaryFileReader
SRCONLY = Main.cpp 
HDRONLY = MixChidist.h My_mvt.h qfc.h Rmath.h WriteLog.h imhof.h Rmath.h

LIB_RAREMETAL_PATH = ../../libRareMetal/
USER_INCLUDES = -I ../../otherLib/eigen-3.2.0 -I$(LIB_RAREMETAL_PATH)include

RMATH_LIB =  $(LIB_RAREMETAL_PATH)/libRmath.a 

# Set additional values prior to including the base Makefiles.
# Set dependencies
opt: $(LIB_RAREMETAL_PATH)libRareMetal.a
debug: $(LIB_RAREMETAL_PATH)libRareMetal_debug.a
profile: $(LIB_RAREMETAL_PATH)libRareMetal_profile.a

# Define USER_LIBS as it should be used in each build type
opt: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal.a $(RMATH_LIB)
debug: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal_debug.a $(RMATH_LIB)
profile: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal_profile.a $(RMATH_LIB)

#OTHER_LIBS = -lRmath

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

USER_COMPILE_VARS = -DVERSION="\"${VERSION}\"" 

binary: USER_COMPILE_VARS += -static
binary: opt
