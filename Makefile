subdir = $(wildcard template/*)
njobs = $(shell cat /proc/cpuinfo | grep "core id" | sort | uniq | wc -l )
expdif = $(wildcard example/*)

default: AFEPack template mpi

AFEPack:
	mkdir -p lib && cd src/ && make -j $(njobs);
	cd ./include && ln -sf ../bundled .;
	@echo "AFEPack done!";

AFEPack_MPI:
	mkdir -p lib && cd mpi/src/ && make;
	@echo "AFEPack MPI done!";

mpi:
	make AFEPack_MPI

template:
	$(foreach dir, $(subdir), cd $(PWD)/$(dir) && make;)
	@echo "template done!";

all: AFEPack template AFEPack_MPI
	@echo "Installation done!";

example:
	$(foreach dir, $(expdif), cd $(PWD)/$(dir) && make;)

clean:
	cd ./src/ && make clean;
	cd ./include/ && rm -rf bundled;
	cd ./mpi/src/ && make clean;
	rm -rf lib;
	$(foreach dir, $(subdir), cd $(PWD)/$(dir) && make clean;)

clean_example:
	$(foreach dir, $(expdif), cd $(PWD)/$(dir) && make clean;)

clean_template:
	$(foreach dir, $(subdir), cd $(PWD)/$(dir) && make clean;)

.PHONY: AFEPack_MPI AFEPack template all example clean clean_example mpi clean_template
