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

using namespace std ;

/* Declaring Global Variables */
GRBEnv env = GRBEnv() ; 
static const int n = 16 ; 
double coeff[n][n] ; 
static const double epsilon = 0.00001; 
static double globalUpperBound = 756980; 
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