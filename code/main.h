#include "gurobi_c++.h"
#include<cassert>
#include<chrono>
#include<queue> 
#include<vector> 
#include<utility>
#include<float.h>
#include<algorithm>
#include<random>
#include<string> 
#include<set> 
#include<fstream>
#include<list>
#include<sys/resource.h>
#include<sys/time.h>
#include<omp.h>
#include<unordered_map>

using namespace std ;

/* Declaring Global Variables */
GRBEnv env = GRBEnv() ; 
static const int n = 20 ; 
static const int max_tabu_size = 5; 
static const int num_iteration = 10 ; 
double coeff[n][n] ; 
static const double epsilon = 0.00001; 
static double globalUpperBound = 1732; 
int a[n][n] ; 
int b[n][n] ; 

/* Necessary Classes */
struct lpProblem
{
    GRBModel model ; 
    double lowerBound ; 
    double upperBound ;
    double* vars ;   

    lpProblem(GRBModel model, double* vars) :
    model(model), lowerBound(0.0), upperBound(0.0), vars(vars) {}

};

/* Function Definitions */

bool double_equal(double, double) ; 
int oneIndex(int, int, int, int) ; 
void manyIndex(int, int*) ; 
lpProblem initialize(const char*) ; 
int lpRelaxation(lpProblem& ) ; 
double findUpperBound(lpProblem&) ;  
void branchAndBound(lpProblem&) ; 
