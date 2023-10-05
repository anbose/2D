CC = mpic++
CFLAGS = -Wall -O3
#Debugging settings
#CFLAGS += -g -DDEBUG
LIBS = -lm
INCLUDE = 

OBJECTS = BEC.o boundary_exchange.o

%.o: %.cc
		${CC} -c ${FLAGS} ${INCLUDES} $< -o $@

all: BEC

BEC: ${OBJECTS}
		mpic++ -o BEC.exe $^ ${LIBS} ${FLAGS}
		
gnuplot:


clean:
		rm -f *.o *.exe

clean-all:
