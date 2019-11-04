include ../../../Makefile.1

CXXFLAGS = -std=c++11 
CXXFLAGS += -I ../../../include/
CXXFLAGS_MPI = $(CXXFLAGS)
CXXFLAGS_MPI += -I ../../../mpi/include/
LDFLAGS = -lboost_serialization
LDFLAGS += -L ../../../lib/ -lAFEPack -ldl
LDFLAGS_MPI = $(LDFLAGS)
LDFLAGS_MPI += -L ../../../lib/ -lAFEPack_mpi

main: main.o data.o SCL.o
	$(MPICXX) -o main main.o data.o SCL.o $(LDFLAGS_MPI)

%.o: %.cpp
	$(MPICXX) -c -o $@ $< $(CXXFLAGS_MPI)

clean:
	-rm -rf main.o data.o SCL.o main
