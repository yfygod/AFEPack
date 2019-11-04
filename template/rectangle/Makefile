include ../../Makefile.1

geometry = $(notdir $(PWD))
subdirs = $(wildcard *.[0-9]) $(wildcard *.[0-9].D)
csourse = $(wildcard *.c)
cppsourse = $(wildcard *.cpp)
object = $(patsubst %.c, %.o, $(csourse)) $(patsubst %.cpp, %.o, $(cppsourse))
library = $(geometry).geometry.so

default: all

all: $(library)
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && cp ../Makefile.tmp Makefile;)
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && $(MAKE) && cp *.bas_fun *.so *.tmp_dof ../;)
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && rm -rf Makefile;)

%.o : %.c
	$(CC) -c -o $@ $< -fPIC -O3

%.o : %.cpp
	$(CXX) -c -o $@ $< -std=c++11 -fPIC -O3

$(library) : $(object)
	$(CC) -shared -o $@ $(object)

clean:
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && cp ../Makefile.tmp Makefile;)
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && make clean;)
	@$(foreach dir, $(subdirs), rm -rf $(dir).bas_fun $(dir).bas_fun.so $(dir).tmp_dof;)
	@$(foreach dir, $(subdirs), cd $(PWD)/$(dir) && rm -rf Makefile;)
	-rm -rf $(object) $(library)
