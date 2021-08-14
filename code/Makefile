PLATFORM = linux64

#Choose correct path when running on HPC 
<<<<<<< HEAD:code/Makefile
TWOUP    = /opt/gurobi912/linux64
=======
TWOUP    = /opt/gurobi/gurobi911/linux64
>>>>>>> main:Makefile
#TWOUP = /opt/apps/gurobi903/gurobi903/linux64

INC      = $(TWOUP)/include/
CC       = gcc
CPP      = g++
<<<<<<< HEAD:code/Makefile
CARGS    = -m64 -g -fopenmp
=======
CARGS    = -m64 -g
>>>>>>> main:Makefile

#Update CLIB and CPPLIB when running on HPC 
CLIB     = -L$(TWOUP)/lib -lgurobi91
CPPLIB   = -L$(TWOUP)/lib -lgurobi_c++ -lgurobi91
#CLIB     = -L$(TWOUP)/lib -lgurobi90
#CPPLIB   = -L$(TWOUP)/lib -lgurobi_c++ -lgurobi90

<<<<<<< HEAD:code/Makefile

all: 
	# Comment out make model when running on HPC 
=======
HELPERDIR = helper

all: 
	# Comment out make model when running on HPC 
	make model 
>>>>>>> main:Makefile
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
