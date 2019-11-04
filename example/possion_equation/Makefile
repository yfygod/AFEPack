source = $(wildcard *.cpp)
maindir = ../../
object = $(patsubst %.cpp, %.o, $(source))
CXXFLAGS = -std=c++11 
INCLUDE += -I $(maindir)/include/
CXXFLAGS += $(INCLUDE) -O3
LDFLAGS = -L $(maindir)/lib -lAFEPack -ldl -pthread


default: main

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(object)
	$(CXX) -o $@ $(object)  $(LDFLAGS) $(LIBS)

clean:
	-rm -rf $(object) main
	-rm -rf u.dx
