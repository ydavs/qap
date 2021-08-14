PLATFORM = linux64

#Choose correct path when running on HPC 
TWOUP    = /opt/gurobi912/linux64
#TWOUP = /opt/apps/gurobi903/gurobi903/linux64

INC      = $(TWOUP)/include/
CC       = gcc
CPP      = g++
CARGS    = -m64 -g -fopenmp

#Update CLIB and CPPLIB when running on HPC 
CLIB     = -L$(TWOUP)/lib -lgurobi91
CPPLIB   = -L$(TWOUP)/lib -lgurobi_c++ -lgurobi91
#CLIB     = -L$(TWOUP)/lib -lgurobi90
#CPPLIB   = -L$(TWOUP)/lib -lgurobi_c++ -lgurobi90


all: 
	# Comment out make model when running on HPC 
	make main

main: main.cpp
	$(CPP) $(CARGS) -g -o main main.cpp -I$(INC) $(CPPLIB) -lm 

cran: cran.cpp
	$(CPP) $(CARGS) -o cran cran.cpp -I$(INC) $(CPPLIB) -lm 

clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name dotnetcore2.csproj -o -name . \) -exec rm -rf {} +; \
	fi
