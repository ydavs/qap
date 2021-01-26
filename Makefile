# Copyright 2020, Gurobi Optimization, LLC

PLATFORM = linux64
TWOUP    = /opt/gurobi/gurobi911/linux64
INC      = $(TWOUP)/include/
CC       = gcc
CPP      = g++
CARGS    = -m64 -g
CLIB     = -L$(TWOUP)/lib -lgurobi91
CPPLIB   = -L$(TWOUP)/lib -lgurobi_c++ -lgurobi91
HELPERDIR = helper

all: 
	make model 
	make main

model: model.o 
	$(HELPERDIR)/genLP

model.o: model.cpp 
	gcc $(HELPERDIR)/genLP.o  -lglpk -lm -o $(HELPERDIR)/genLP

model.cpp:
	gcc -c $(HELPERDIR)/genLP.c  -o $(HELPERDIR)/genLP.o

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