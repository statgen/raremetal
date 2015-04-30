VERSION=4.13
EXE ?= libRareMetal.a

CURRENT_PATH := $(dir $(lastword $(MAKEFILE_LIST)))
BINDIR=$(CURRENT_PATH)
BINDIR_DEBUG= 
BINDIR_PROFILE= 
PROG_OPT = $(BINDIR)$(EXE)
PROG_DEBUG = $(PROG_OPT:.a=_debug.a)
PROG_PROFILE = $(PROG_OPT:.a=_profile.a)

OBJDIR = obj

########################
# Library Path:
#
# Default path to the statgen library.  
# If you want to use a different path, you can:
#   1) change this file
#   2) set LIB_PATH_GENERAL to the location all libStatGen tools should use.
#   3) set LIB_PATH_LIB_RAREMETAL to the location libRareMetal should use.
LIB_PATH_GENERAL ?=../libStatGen
LIB_PATH_LIB_RAREMETAL ?= $(LIB_PATH_GENERAL)

# add any additional ../ as necessary if it is a relative path
ACTUAL_PATH := $(patsubst ../%, $(CURRENT_PATH)../%, $(LIB_PATH_LIB_RAREMETAL))
INCLUDE_MAKE := $(ACTUAL_PATH)/Makefiles/$(PARENT_MAKE)

########################
# Include the base Makefile
include $(INCLUDE_MAKE)

.phony: cloneLib

########################
# This is to handle the case where the statgen library is not already there.
$(INCLUDE_MAKE):
	@echo Unable to locate: $(INCLUDE_MAKE)
	@echo To change the location, set LIB_PATH_GENERAL or LIB_PATH_LIB_RAREMETAL to the appropriate path to libStatGen.  Or specify \"make LIB_PATH_GENERAL=yourPath\" or \"make LIB_PATH_LIB_RAREMETAL=yourPath\"
	@echo Use make cloneLib if you have git and want to clone the current libStatGen at that location.

cloneLib:
	@if test -d $(LIB_PATH_LIB_RAREMETAL); \
	then echo $(LIB_PATH_LIB_RAREMETAL) already exists; \
	else git clone git://github.com/statgen/libStatGen.git $(LIB_PATH_LIB_RAREMETAL); fi
	@echo Call make to compile libStatGen and this tool.

########################
# add headers into include directory.
EXTLIB_INCLUDE_PATH=$(CURRENT_PATH)include

# make everything, ensure headers are in the include direcotry.
opt debug profile : $(addprefix $(EXTLIB_INCLUDE_PATH)/, $(TOOLHDR)) 

# Link into the include directory
$(EXTLIB_INCLUDE_PATH)/%.h: %.h
	ln  -s ../$(DIR_NAME)/$< $@

USER_INCLUDES = -I$(EXTLIB_INCLUDE_PATH)
USER_REMOVES = -rm -f $(addprefix $(EXTLIB_INCLUDE_PATH)/, $(TOOLHDR)) 
