include ../../Makefile.1

CXXFLAGS = -std=c++11 
CXXFLAGS += -I ../../include/
CXXFLAGS_MPI = $(CXXFLAGS)
CXXFLAGS_MPI += -I ../../mpi/include/
LDFLAGS = -lboost_serialization
LDFLAGS += -L ../../lib/ -lAFEPack -ldl
LDFLAGS_MPI = $(LDFLAGS)
LDFLAGS_MPI += -L ../../lib/ -lAFEPack_mpi
bin = $(patsubst %.cpp, %, $(wildcard *.cpp)) 

default: $(bin) ex19 ex20 

ex03: ex03.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(LDFLAGS)

ex04: ex04.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(LDFLAGS)

ex05: ex05.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex06: ex06.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex07: ex07.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex08: ex08.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex09: ex09.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex10: ex10.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex11: ex11.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex12: ex12.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex13: ex13.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex14: ex14.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex15: ex15.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex16: ex16.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex18: ex18.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex19: 
	cd ex19 && make

ex20: 
	cd ex20 && make

ex28: ex28.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex29-2d: ex29-2d.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

ex29-3d: ex29-3d.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)


ex30: ex30.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS_MPI) $(LDFLAGS_MPI)

.PHONY: ex19 ex20

clean:
	-rm -rf $(bin)
	cd ex19 && make clean
	cd ex20 && make clean
