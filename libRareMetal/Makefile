SUBDIRS=mvtnorm pdf src

PARENT_MAKE := Makefile.tool
include Makefile.inc

USER_LIBS = libRmath.a

pdf: src

clean:$(SUBDIRS)
	rm -f $(PROG_OPT)
	rm -f $(PROG_DEBUG)
	rm -f $(PROG_PROFILE)
