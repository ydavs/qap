/* QAP  */

#include "main.h"

/* Helper Functions */
bool double_equal(double a, double b)
{
  return abs(a-b) < epsilon ; 
}

int oneIndex(int i, int j, int k, int l)
{
    return (i-1)*n*n*n + (j-1)*n*n + (k-1)*n + (l-1) ; 
}

void manyIndex(int val, int* array)
{
    array[0] = (int) val/(n*n*n) + 1 ; 
    val = val - (array[0]-1)*n*n*n ; 

    array[1] = (int) val/(n*n) + 1 ; 
    val = val - (array[1] - 1)*n*n ; 

    array[2] = (int) val/n + 1 ; 
    val = val - (array[2] - 1)*n ; 

    array[3] = (int) val + 1 ;
    val = val - (array[3] -1) ; 

    //cout << val << " " << array[0] << " " << array[1] << " " << array[2] << " " << array[3] << endl; 
    assert(val == 0) ; 

    return ;
}

void readCoeff(char* fname)
{
    ifstream file ; 
    file.open(fname) ; 

    double dump ; 
    file >> dump ; 

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)file >> a[i][j] ; 
    
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++) file >> b[i][j] ; 
    
    return; 
    
}

/* Helper Functions */


/* Initializations */
lpProblem initialize(const char* fname)
{
    GRBModel model = GRBModel(env, fname) ; 
    static double arr[n*n*n*n] ; 

    /* Setting Model Parameters */
    model.set(GRB_IntParam_OutputFlag, 1) ; 
    model.set(GRB_StringParam_LogFile, "logfile.log") ; 
    model.set(GRB_IntParam_LogToConsole, 0) ; 

    lpProblem prob(model, arr) ; 

    return prob ; 
}

lpProblem initialize2(GRBModel& model)
{
    static double arr[n*n*n*n] ; 

    /* Setting model parameters */
    model.set(GRB_IntParam_OutputFlag, 1) ; 
    model.set(GRB_StringParam_LogFile, "logfile.log") ; 
    model.set(GRB_IntParam_LogToConsole, 0) ; 
    
    lpProblem prob(model, arr) ; 

    return prob ; 
}

/* Initializations */

/* QUBO Solver */
int qubo_solver(vector<int> &tmp)
{
  try
  {
    GRBModel model = GRBModel(env) ; 
    GRBQuadExpr objExp ; 
    GRBLinExpr cons ; 
    double rngcoeff[n] ;
    GRBVar var[n] ;   

    //cout << "hel" << endl; 
    model.set(GRB_IntParam_OutputFlag, 0) ; 

    /*  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++) cout << coeff[i][j] << " " ; 
      cout << endl ; 
    }  */

    for(int i=0;i<n;i++)
    {
      var[i] = model.addVar(0.0, 1.0, GRB_INFINITY, GRB_BINARY, to_string(11)) ;
      rngcoeff[i] = 1 ;  
    }
    
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
      {
        if(i >= j) continue ; 
        objExp.addTerm(coeff[i][j], var[i], var[j]) ; 
      }
    }

    model.setObjective(objExp, GRB_MAXIMIZE) ; 
    
    model.optimize() ; 

    int optimstatus = model.get(GRB_IntAttr_Status) ; 

    if(optimstatus == GRB_INF_OR_UNBD)
    {
      model.set(GRB_IntParam_Presolve, 0);
      model.optimize();
      optimstatus = model.get(GRB_IntAttr_Status);
    }

    if (optimstatus == GRB_OPTIMAL) 
    {
      //double objval = model.get(GRB_DoubleAttr_ObjVal);
      //cout << "Optimal objective: " << objval << endl;

      for(int i=0;i<n;i++)
      {
        if(var[i].get(GRB_DoubleAttr_X)) tmp.push_back(i+1) ; 
      }
    }
    else if (optimstatus == GRB_INFEASIBLE) 
    {
      cout << "Model is infeasible" << endl;

      // compute and write out IIS
      model.computeIIS();
      model.write("model.ilp");

      return -1 ; 
    }
    else if (optimstatus == GRB_UNBOUNDED) 
    {
      cout << "Model is unbounded" << endl;
    } 
    else 
    {
      cout << "Optimization was stopped with status = "
           << optimstatus << endl;
    }
  }
  catch(GRBException e) 
  {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch (...) 
  {
    cout << "Error during optimization" << endl;
  }

  return 0 ; 
}

/* Linear Relaxation */
int lpRelaxation(lpProblem& prob)
{
    try
    {
        prob.model.set(GRB_DoubleParam_IterationLimit, 30000) ; 
        prob.model.set(GRB_IntParam_Method, 2) ; 
        prob.model.optimize() ; 

        int optimstatus = prob.model.get(GRB_IntAttr_Status) ; 
        if(optimstatus == GRB_INF_OR_UNBD)
        {
            prob.model.set(GRB_IntParam_Presolve, 0) ; 
            prob.model.optimize() ; 
            optimstatus = prob.model.get(GRB_IntAttr_Status) ; 
        }

        if(optimstatus == GRB_OPTIMAL)
        {
            double objval = prob.model.get(GRB_DoubleAttr_ObjVal) ; 
            prob.lowerBound = objval ; 

            /* Storing values of variable for optimal value */
            GRBVar* vars = prob.model.getVars() ; 
            for(int i=0;i<n*n*n*n;i++) prob.vars[i] = vars[i].get(GRB_DoubleAttr_X) ;
        } else if(optimstatus == GRB_INFEASIBLE)
        {
            cout << "Model is infeasible" << endl ; 
            /* Compute and write out IIS */
            prob.model.computeIIS() ; 
            prob.model.write("model.ilp") ; 
            return 1 ; 
        } else if(optimstatus == GRB_UNBOUNDED){ cout << "Model is unbounded " << endl ; return 1 ; }
        else cout << "Optmization is stopped with status" << optimstatus << endl ;
    }
    catch(GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl ; 
        cout << e.getMessage() << endl ; 
    } catch(...)
    {
        cout << "Error during optmization" << endl ; 
    }
    return 0 ; 
}

/* Upper Bound Calculation  */
double findUpperBound(lpProblem& prob)
{
    bool Itaken[n+1], Jtaken[n+1] ; 
    int myJ[n+1] ; 
    double val = 0 ; 

    for(int i=0;i<=n;i++) Itaken[i] =0, Jtaken[i] =0, myJ[i] =0 ; 

    vector<pair<double, pair<int, int>>> lst ; 

    for(int i=1;i<=n;i++)
    {
        for(int j=1;j<=n;j++) lst.push_back({prob.vars[oneIndex(i, j, i ,j)], {i,j}}) ; 
    }

    sort(lst.begin(), lst.end(), greater<>()) ; 

    for(auto i : lst)
    {
        if(!Itaken[i.second.first] && !Jtaken[i.second.second])
        {
            myJ[i.second.first] = i.second.second ; 
            Itaken[i.second.first] = 1 ; 
            Jtaken[i.second.second] = 1 ; 
        }
    }

    /* for(int i : myJ) cout << i << " " ; 
    cout << endl ;  */

    for(int i=1;i<=n;i++)
    {
        for(int j=1;j<=n;j++) val+= a[i-1][j-1]*b[myJ[i]-1][myJ[j]-1] ; 
    }

    prob.upperBound = val ; 
    return prob.upperBound ; 
}

/* Cutting Procedures */
/* Kaibel Box 1 Cuts */
void box1cuts(lpProblem& prob)
{
    random_device rd ; 
    mt19937 gen(rd()) ; 
    uniform_int_distribution<> distrib(1, n) ;

    GRBVar* vars = prob.model.getVars() ; 
    int cuts = 0 ; 

    for(int it=1;it<=2;it++)
    {
        cout << "Iteration : " << it << endl ; 
        int cuts = 0 ; 
        int cntlimit = 1000 ; 

        for(int beta=2;beta<4;beta++)
        {
            for(int ix=0;ix<cntlimit;ix++)
            {
                int qsize = n/2 +2 -it ;
                //int qsize = qdis(gen) ;  
                set<int> q ;

                for(int i=0;i<qsize;i++)
                {
                    int ins = distrib(gen) ;
                    while(q.find(ins)!=q.end() && q.size()) ins = distrib(gen) ; 
                    q.insert(ins) ;
                }

                for(int a=0;a<n;a++)
                {
                    double fifth = 0.0 ; 
                    for(int i : q) fifth += (beta-1)*prob.vars[oneIndex(a+1, i, a+1, i)] ; 

                    for(int k=0;k<n;k++)
                    {
                        double first = 0.0 ; 

                        if(a==k){coeff[a][k] =  -10 ; continue ;}

                        for(int i :q)
                        {
                            for(int j : q) first+= prob.vars[oneIndex(a+1,i,k+1,j)] ; 
                        }
                        coeff[a][k] = -0.8 + fifth - first ; 
                    }
                }

                vector<int> p ; 
                if(qubo_solver(p) == -1) continue ; 

                assert(qsize >= beta+1 && qsize <= n-3) ; 
                assert(qsize+p.size() <= n-3+beta) ; 

                GRBLinExpr exp  ;
                double val =0.0, pos = 1.0, neg = -1.0 ; 

                for(int i : p)
                {
                    for(int j :q)
                    {
                        val+=(beta-1)*prob.vars[oneIndex(i,j,i,j)] ; 
                        double tmpVal = beta-1;
                        exp.addTerms(&tmpVal, &vars[oneIndex(i,j,i,j)], 1) ; 

                        for(int k : p)
                        {
                            if(i>=k) continue ; 

                            for(int l : q)
                            {
                                val-= prob.vars[oneIndex(i,j,k,l)] ; 
                                exp.addTerms(&neg, &vars[oneIndex(i,j,k,l)], 1) ; 
                            }
                        } 
                    }
                }

                double limit = (double) (beta*beta - beta)/2 ; 
                if(val > limit + 0.05)
                {
                    prob.model.addConstr(exp, GRB_LESS_EQUAL, limit) ; 
                    cuts++ ; 
                }
            }

            for(int ix=0;ix<cntlimit;ix++)
            {
                int qsize = n/2 + 2 - it; 
                //int qsize = qdis(gen) ; 
                set<int> q ;

                for(int i=0;i<qsize;i++)
                {
                    int ins = distrib(gen) ;
                    while(q.find(ins)!=q.end() && q.size()) ins = distrib(gen) ; 
                    q.insert(ins) ;
                }

                for(int a=0;a<n;a++)
                {
                    double fifth = 0.0 ; 
                    for(int i : q) fifth += (beta-1)*prob.vars[oneIndex(i, a+1, i, a+1)] ; 

                    for(int k=0;k<n;k++)
                    {
                        double first = 0.0 ; 

                        if(a==k){coeff[a][k] =  -10 ; continue ;}

                        for(int i :q)
                        {
                            for(int j : q) first+= prob.vars[oneIndex(i,a+1,j,k+1)] ; 
                        }
                        coeff[a][k] = -0.8 + fifth - first ; 
                    }
                }

                vector<int> p ; 
                if(qubo_solver(p) == -1) continue ; 

                assert(qsize >= beta+1 && qsize <= n-3) ; 
                assert(qsize+p.size() <= n-3+beta) ; 

                GRBLinExpr exp  ;
                double val =0.0, pos = 1.0, neg = -1.0 ; 

                for(int i : p)
                {
                    for(int j :q)
                    {
                        val+=(beta-1)*prob.vars[oneIndex(j,i,j,i)] ; 
                        double tmpVal = beta-1;
                        exp.addTerms(&tmpVal, &vars[oneIndex(j,i,j,i)], 1) ; 

                        for(int k : p)
                        {
                            if(i>=k) continue ; 

                            for(int l : q)
                            {
                                val-= prob.vars[oneIndex(j,i,l,k)] ; 
                                exp.addTerms(&neg, &vars[oneIndex(j,i,l,k)], 1) ; 
                            }
                        } 
                    }
                }

                double limit = (double) (beta*beta - beta)/2 ; 
                if(val > limit + 0.05)
                {
                    prob.model.addConstr(exp, GRB_LESS_EQUAL, limit) ; 
                    cuts++ ; 
                }
            }
        }

        cout << "Cuts added : " << cuts << endl ; 
    }
}

/* QAP2 Cuts */
void qap2(lpProblem& prob)
{
    int m = 3 ;
    int cuts = 0 ; 
    GRBVar* vars = prob.model.getVars() ; 

    for(int m=3;m<=n;m++){
    for(int k=1;k<=n;k++)
    {
        for(int l=1;l<=n;l++)
        {
            //cout << k << " " << l << " " << endl ; 
            vector<pair<double, pair<int, int>>> lst ; 
            
            for(int i=1;i<=n;i++)
                for(int j=1;j<=n;j++) lst.push_back({prob.vars[oneIndex(i,j,k,l)], {i,j}}) ; 

            sort(lst.begin(), lst.end(), greater<>()) ;

            bool itaken[n+1], jtaken[n+1] ; 

            for(int i=0;i<=n;i++) itaken[i] = 0, jtaken[i] = 0 ; 
            itaken[k] = 1, jtaken[l] = 1 ; 

            vector<int> ivec, jvec ; 

            for(auto i : lst)
            {
                if(ivec.size() > m-1 && jvec.size() > m-1) break ; 
                if(!itaken[i.second.first] && !jtaken[i.second.second])
                {
                    ivec.push_back(i.second.first) ; 
                    jvec.push_back(i.second.second) ; 
                }
            }  

            GRBLinExpr exp ; 
            double first = 0.0, third = 0.0 ; 
            double tmp1 =1.0, tmp2 = -1.0 ; 
            for(int i=0;i<m;i++)
            { 
                first += prob.vars[oneIndex(ivec[i], jvec[i], k,l)] ;
                exp.addTerms(&tmp1, &vars[oneIndex(ivec[i], jvec[i], k, l)], 1) ;  
                for(int j=0;j<i;j++){ 
                    third += prob.vars[oneIndex(ivec[j], jvec[j], ivec[i], jvec[i])] ; 
                    exp.addTerms(&tmp2, &vars[oneIndex(ivec[j], jvec[j], ivec[i], jvec[i])], 1) ; 
                }
            }

            double sum = first - prob.vars[oneIndex(k, l, k,l)] - third ; 
            exp.addTerms(&tmp2, &vars[oneIndex(k,l,k,l)], 1) ; 

            if(sum > 0.05)
            {
                prob.model.addConstr(exp, GRB_LESS_EQUAL, 0.0) ;
                //cout << sum << endl ; 
                cuts++ ;  
            }
        }
    } }

    cout << "Cuts Added : " << cuts << endl; 
    return  ; 
}


/* Branch and bound method */
void branchAndBound(lpProblem& prob)
{
    queue<lpProblem> q  ; 
    q.push(prob) ; 

    int it = 0 ;  
    while(!q.empty())
    {
        if(it >0 ) break ; 
        it++ ; 
        cout << it << " Global Lower Bound : " << globalLowerBound << " Global Upper Bound : " << globalUpperBound << endl ;  

        if(lpRelaxation(q.front())) {q.pop(); continue;} 
        findUpperBound(q.front()) ; 

        cout << q.front().lowerBound << " " << q.front().upperBound << endl ; 

        globalUpperBound = min(globalUpperBound, q.front().upperBound)  ; 
        if(double_equal(q.front().lowerBound, q.front().upperBound) || q.front().lowerBound >= globalUpperBound)
        {
            //globalLowerBound = max(q.front().lowerBound, globalLowerBound) ; 
            q.pop() ; 
            continue ; 
        }

        /*  Cutting Procedures  */
        qap2(q.front()) ;
        //box1cuts(q.front()) ;  

        /* Resolving after adding Cuts */
        lpRelaxation(q.front()) ; 
        findUpperBound(q.front()) ; 

        globalUpperBound = min(globalUpperBound, q.front().upperBound)  ; 
        if(q.front().lowerBound >= globalUpperBound){ q.pop() ; continue ; } 
        if(double_equal(q.front().lowerBound, q.front().upperBound))
        {
            //globalLowerBound = max(q.front().lowerBound, globalLowerBound) ; 
            q.pop() ; 
            continue ; 
        }

        cout << "Adding Box1 Cuts(Trivial Implementation)" << q.front().lowerBound << " " << q.front().upperBound << endl ; 


        /* Branching on max value found on diagonal */
        double val = 0 ; 
        pair<int, int> p ;  
        for(int i=1;i<=n;i++)
        {
            for(int j=1;j<=n;j++)
            {
                if(q.front().vars[oneIndex(i,j,i,j)] > val && !double_equal(q.front().vars[oneIndex(i,j,i,j)], 1))
                {
                    val = q.front().vars[oneIndex(i,j,i,j)] ; 
                    p = {i,j} ; 
                }
            }
        }

        cout << "Node for branching is : " << p.first << " " << p.second  << " " << val << endl ; 

        /* Adding branching node's binary values to the model */
        for(int i=1;i>=0;i--)
        {
            GRBModel model = GRBModel(q.front().model) ; 
            lpProblem newProb = initialize2(model) ; 

            if(i) newProb.model.getVar(oneIndex(p.first, p.second, p.first, p.second)).set(GRB_DoubleAttr_LB, 1.0) ;
            else newProb.model.getVar(oneIndex(p.first, p.second, p.first, p.second)).set(GRB_DoubleAttr_UB, 0.0) ; 

            newProb.model.update() ; 

            q.push(newProb) ; 
        }
        q.pop() ; 
    }
    return ; 
}

int main(int argc, char* argv[])
{
    auto start = chrono::steady_clock::now() ; 
    cout.precision(17) ; 

    if(argc < 3)
    {
        cout << "Usage : ./main model instance" << endl ; 
        return 1;  
    }

    readCoeff(argv[2]) ; 
    lpProblem prob = initialize(argv[1]) ; 
    branchAndBound(prob) ; 

    auto end = chrono::steady_clock::now() ; 
    cout << "Time Taken :  " << (double) chrono::duration<double, milli>(end - start).count()/60000 << "min" << endl ; 

    return 0 ; 
}