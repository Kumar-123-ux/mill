# Project: mill
# Makefile created by Dev-C++ 4.9.9.2

gcc  = gcc 
g++   = g++
WINDRES = windres.exe
RES  = 
OBJ  =  dem.o init.o injection.o facesort.o fileio.o motion.o  contact.o force.o mathop.o neighbour.o$(RES)
LINKOBJ  = dem.o init.o injection.o facesort.o fileio.o motion.o  contact.o force.o mathop.o neighbour.o$(RES)
LIBS =  -lm -L"lib"   
INCS =  -I"include" 
CXXINCS =  -I"include" 
BIN  = mill 
CXXFLAGS = $(CXXINCS) -g3  
CFLAGS = -O -systype bsd43 $(INCS)  -g3
RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all-before:
all-after:
clean-custom:

all: all-before mill all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(gcc) $(LINKOBJ) -o "mill"  $(LIBS)

dem.o: dem.c
	$(gcc) -c dem.c -o dem.o $(CXXFLAGS)

init.o: init.c
	$(gcc) -c init.c -o init.o $(CXXFLAGS)

injection.o: injection.c
	$(gcc) -c injection.c -o injection.o $(CXXFLAGS)

facesort.o: facesort.c
	$(gcc) -c facesort.c -o facesort.o $(CXXFLAGS)

fileio.o: fileio.c 
	$(gcc) -c fileio.c -o fileio.o $(CXXFLAGS)

motion.o: motion.c
	$(gcc) -c motion.c -o motion.o $(CXXFLAGS)

contact.o: contact.c
	$(gcc) -c contact.c -o contact.o $(CXXFLAGS)
force.o: force.c
	$(gcc) -c force.c -o force.o $(CXXFLAGS)
mathop.o: mathop.c
	$(gcc) -c mathop.c -o mathop.o $(CXXFLAGS)
neighbour.o: neighbour.c
	$(gcc) -c neighbour.c -o neighbour.o $(CXXFLAGS)
