PROGRAM = gen
OBJS	= main.o
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
