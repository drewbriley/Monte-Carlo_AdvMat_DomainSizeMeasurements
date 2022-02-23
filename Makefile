
rho ?= 0
all: programs

INCLUDE=./include
%.o: ./src/%.cpp Makefile
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -O3 -fopenmp -c $< -o $@ -lstdc++fs

# Programs section
programs: 3D-MC-DomainSizeExp_$(rho)

3D-MC-DomainSizeExp_$(rho): tools.o \
          Exciton.o \
          Lattice_3D.o \
          MonteCarloExperiments.o \
          3D-MC-DomainSizeExp.o
	$(CXX) $(CXXFLAGS) -fopenmp $^ -o $@ -lstdc++fs


