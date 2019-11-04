source = $(wildcard *.cpp)
maindir = ../../
object = $(patsubst %.cpp, %.o, $(source))
exebin = $(patsubst %.cpp, %, $(source))
CXXFLAGS = -std=c++11 
INCLUDE += -I $(maindir)/include/
CXXFLAGS += $(INCLUDE)
LDFLAGS = -L $(maindir)/lib -lAFEPack -ldl -pthread


default: main

order=2

% : %.cpp
	$(CXX)  -o $@ $< $(CXXFLAGS) $(LDFLAGS)

main: $(exebin)

clean:
	-rm  $(exebin)
