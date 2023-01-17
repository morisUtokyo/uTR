PROGRAM = uTR
OBJS	= main.o handle_one_file.o SAIS.o coverage_by_units.o units.o coverage_by_long_units_nsop_Z.o string_decomposer.o smooth.o
CC	= gcc
CPP	= g++
CFLAGS	= -std=c99 -fPIC -fcommon

.cpp.o:
	$(CPP) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<

# g++ must be used to link libraries required
$(PROGRAM): $(OBJS)
	$(CPP) $(OBJS) -o $(PROGRAM)
clean:
	rm $(PROGRAM) $(OBJS)
