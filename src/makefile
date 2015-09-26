#
# Makefile for Simple MC
# John R. Dowdle
# June 2008
#
#
#
CC = mpicxx
LD = mpicxx
# CFLAGS = -Wall -pedantic -O3 -Wno-long-long -g -pg
# LDFLAGS = -lxdrfile -lgsl -lgslcblas -g -pg
CFLAGS = -Wall -pedantic -O3 -Wno-long-long
LDFLAGS = -lxdrfile -lgsl -lgslcblas

PROG = jls
OBJS = main.o random_seed.o particle_system_construct.o particle_system_mc.o particle_system_free_eng.o particle_system_fileio.o particle_system_utils.o particle_system_init.o particle_system_pot_eng.o utils.o particle.o multi_system.o cell.o cell_list.o

all:	$(PROG)

$(PROG):	$(OBJS)
		$(LD) $(OBJS) $(LDFLAGS) -o $(PROG)

particle_system_%.o:	particle_system_%.cpp particle_system.h
	$(CC) $(CFLAGS) -c $<

%.o:	%.cpp %.h
	$(CC) $(CFLAGS) -c $<

main.o:		main.cpp random_seed.h particle_system.h utils.h particle.h
	$(CC) $(CFLAGS) -c $<

clean:
	/bin/rm -f *.o jls

clobber:
	/bin/rm -rf *~ \#*\# *.o $(PROG) test* core

