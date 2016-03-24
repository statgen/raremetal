VERSION=4.13.8

CURRENT_PATH := $(dir $(lastword $(MAKEFILE_LIST)))
BINDIR=$(CURRENT_PATH)bin
LIB_PATH_GENERAL ?=libStatGen
LIB_PATH_RMW ?= $(LIB_PATH_GENERAL)

# add any additional ../ as necessary if it is a relative path
ifeq ($(LIB_PATH_RMW:/%=/),/)
   ACTUAL_PATH = $(LIB_PATH_RMW)
else
   ACTUAL_PATH = $(patsubst %, $(CURRENT_PATH)%, $(LIB_PATH_RMW))
endif
INCLUDE_MAKE = $(ACTUAL_PATH)/Makefiles/$(PARENT_MAKE)

########################
# Include the base Makefile
include $(INCLUDE_MAKE)

.phony: cloneLib

$(INCLUDE_MAKE):
	@echo Unable to locate: $(INCLUDE_MAKE)
	@echo To change the location, set LIB_PATH_GENERAL or LIB_PATH_RMW to the appropriate path to libStatGen.  Or specify \"make LIB_PATH_GENERAL=yourPath\" or \"make LIB_PATH_RMW=yourPath\"
	@echo Use make cloneLib if you have git and want to clone the current libStatGen at that location.

cloneLib:
	@if test -d $(LIB_PATH_RMW); \
	then echo $(LIB_PATH_RMW) already exists; \
	else git clone git://github.com/statgen/libStatGen.git $(LIB_PATH_RMW); fi
	@echo Call make to compile libStatGen and this tool.

WHOLEPACKAGE_CMD = tar chvz --exclude="*~" --exclude=$(RELEASE_FILE) --exclude='obj/*' --exclude='libRareMetal.a' --exclude='libRareMetal_debug.a' --exclude='libRareMetal_profile.a' --exclude='libStatGen.a'  --exclude='libStatGen_debug.a'  --exclude='libStatGen_profile.a' --exclude='include/*' --exclude='bin/*' --exclude='test/results/*' --exclude-vcs --transform 's,^,$(DIR_NAME)_$(VERSION)/,' -C .. $(DIR_NAME) -C $(DIR_NAME) -C $(DIR_ABOVE_LIB) $(BASE_LIB_DIRNAME) --show-transformed-names

wholepackageLinux:
	touch LINUX_$(RELEASE_FILE)
	$(WHOLEPACKAGE_CMD) -f LINUX_$(RELEASE_FILE) --exclude=*_$(RELEASE_FILE)

wholepackageMac:
	touch MAC_$(RELEASE_FILE)
	$(WHOLEPACKAGE_CMD) -f MAC_$(RELEASE_FILE) --exclude=*_$(RELEASE_FILE) --transform 's,libRmath_MAC.a,libRmath.a,'

wholepackageCygwin64:
	touch CYGWIN64_$(RELEASE_FILE)
	$(WHOLEPACKAGE_CMD) -f CYGWIN64_$(RELEASE_FILE) --exclude=*_$(RELEASE_FILE) --transform 's,libRmath_cygwin64.a,libRmath.a,' 

wholepackageMingw:
	touch MINGW_$(RELEASE_FILE)
	$(WHOLEPACKAGE_CMD) -f MINGW_$(RELEASE_FILE) --exclude=*_$(RELEASE_FILE) --transform 's,libRmath_mingw.a,libRmath.a,'

wholepackageAll: wholepackageLinux wholepackageMac wholepackageCygwin64 wholepackageMingw
	@echo "Made all packages"


