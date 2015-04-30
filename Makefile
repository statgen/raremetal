SUBDIRS = raremetal raremetalworker libRareMetal 
raremetal raremetalworker: libRareMetal

PARENT_MAKE := Makefile.tool
include Makefile.inc

USER_WHOLEPACKAGE_DIRS = libRareMetal
