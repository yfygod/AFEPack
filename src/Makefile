include ../Makefile.1

indir = -I ../include/
source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
CXXFLAGS += -std=c++11 -O3 -march=native -D__SERIALIZATION__ -fPIC -g 
#CXXFLAGS += -DENABLE_AFEPACK_NAMESPACE
all : libAFEPack.so

libAFEPack.so : $(object)
	$(CXX) -shared -o ../lib/$@ $(object)

%.o : %.cpp
	$(CXX) -c -o $@ $< $(indir) $(CXXFLAGS)

clean :
	-rm -rf $(object) ../lib/libAFEPack.so

.PHONY : default clean all


