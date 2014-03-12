UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
FC = gfortran
endif
ifeq ($(UNAME), Darwin)
FC = gfortran
endif
#FC = ~/Downloads/gcc-4.6/bin/gfortran
DEBUGFLAGS = -g3 -pg -O0 \
	-fbacktrace -fbounds-check -Wall
RUNFLAGS = -O3
OBJDIRBASE := ./obj/
BINDIR = ./bin/
LIBDIR = ./lib/
DATDIR = ./dat/
SRCDIR = ./src/
TESTDIR= $(SRCDIR)Test/
INPUTDIR=./input/
DEBUGDIR=$(OBJDIRBASE)debug/
RUNDIR = $(OBJDIRBASE)run/
FFLAGS = -I $(OBJDIR) 
LIBS = recipes_f90
LDFLAGS= -l $(LIBDIR)$(LIBS)
vpath %.f90 $(SRCDIR):$(TESTDIR)
vpath %.o $(OBJDIR)
vpath %.mod $(OBJDIR)
%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
	mv $@ $(OBJDIR)
	if test -e $(@:o=mod) ; then mv $(@:o=mod) $(OBJDIR); fi

OBJS = $(addprefix $(OBJDIR),$(addsuffix .o, global_data types tester))


.PHONY : clean purge debug run

debug : FFLAGS += $(DEBUGFLAGS) 
debug : OBJDIR = $(OBJDIRBASE)debug/
debug : build testbin
build : $(OBJS)
testbin : tester.o build
	$(FC) $(FFLAGS)-o $@ $(OBJS)
	mv $@ $(BINDIR)

clean :
	-rm -r $(BINDIR)*.* $(DEBUGDIR)*.* $(RUNDIR)*.*
purge : clean
	-rm -r $(DATDIR)*.*

types.o  : global_data.o
tester.o : types.o
