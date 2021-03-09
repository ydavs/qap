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

            /* Removing all loose constraints  */
            GRBConstr* constrs = prob.model.getConstrs() ; 
            //cout << "Constraint Number bf : " << prob.model.get(GRB_IntAttr_NumConstrs) << endl; 

            for(int i=0;i<prob.model.get(GRB_IntAttr_NumConstrs);i++)
                if(constrs[i].get(GRB_DoubleAttr_Slack) != 0 ) prob.model.remove(constrs[i]) ;
            prob.model.update() ; 
            
            //cout << "Constraint Number af : " << prob.model.get(GRB_IntAttr_NumConstrs) << endl; 
            delete vars, constrs ; 
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

/* tabu Search */
double fitness(lpProblem &prob, vector<int> &p, vector<int> &q, int cut, int beta)
{
     /* 
        Cut = 0 -> tabu search for Box1 inequalities
        Cut = 1 -> tabu search for QAP2 inequalities
     */
    if(!cut)
    {
        GRBVar* vars = prob.model.getVars() ; 
        double val=0.0 ; 

        for(int i : p)
        {
            for(int j :q)
            {
                val+=(beta-1)*prob.vars[oneIndex(i, j,i,j)] ; 

                for(int k : p)
                {
                    if(i>=k) continue ; 
                    for(int l :q)
                    {
                        val-=prob.vars[oneIndex(i,j,k,l)] ; 
                    }
                }
            }
        }
        delete vars ; 
        return val ;
    }
    else
    {
        GRBVar* vars = prob.model.getVars() ; 
        double val=0.0 ; 

        assert(q.size() == p.size()) ; 
        int k = p[p.size()-1] ; 
        int l = q[p.size()-1] ; 

        val -= prob.vars[oneIndex(k,l,k,l)] ; 
        for(int i=0;i<p.size()-1;i++)
        {
            val+= prob.vars[oneIndex(p[i], q[i], k,l)] ; 
            for(int j=0;j<p.size()-1;j++)
            {
                if(j<=i) continue ; 
                val-= prob.vars[oneIndex(p[i], q[i], p[j], q[j])] ; 
            }
        }

        delete vars ; 
        return val ; 
    } 
}

void get_neighbours(list<vector<int>> &que, vector<int> &p)
{
    bool bit_array[n+1] ; 
    for(int i=0;i<=n;i++) bit_array[i] = 0 ; 
    for(int i : p ) bit_array[i] = 1;  
    for(int i=0;i<p.size();i++)
    {
        vector<int> tmp(p.size()) ; 
        for(int j=0;j<p.size();j++)
        {
            if(j == i) continue ;
            tmp[j] = p[j] ; 
        }

        for(int j=1;j<=n;j++)
        {
            if(!bit_array[j]) tmp[i] = j ; 
            else continue ;
            que.push_back(tmp) ; 
        }
    }
}

bool check_list(list<vector<int>> &lst, vector<int> &vec)
{
    for(vector<int> candidate : lst)
    {
        assert(candidate.size() == vec.size()) ; 
        bool chk = 1 ; 
        for(int i=0;i<vec.size();i++)
        {
            if(candidate[i] != vec[i]) chk = 0 ; 
        }
        if(chk==1) return 1 ; 
    }
    return 0 ; 
}

void tabu_search(lpProblem &prob, vector<int> &p, vector<int> &q,  int psize, int cut, int beta)
{  
    vector<int> best = p ; 
    vector<int> best_candidate = p ;
    list<vector<int>> tabulist ; 
    tabulist.push_back(p) ; 
    int it = 0, max_tabu_size=20, num_iteration=40 ; 

    //if(cut) num_iteration = 5 ; 

    while(it < num_iteration) 
    {
        it++ ; 
        list<vector<int>> neighbours ;  
        get_neighbours(neighbours, best_candidate) ; 
        best_candidate = neighbours.front() ;

        for(vector<int> candidate : neighbours)
        {
            if(!check_list(tabulist, candidate) && (fitness(prob,candidate, q, cut, beta) > fitness(prob, best_candidate, q, cut, beta))) best_candidate = candidate ;
            if(fitness(prob, best_candidate, q, cut, beta) > fitness(prob, best, q, cut, beta)) best = best_candidate ; 
        }
        tabulist.push_back(best_candidate) ; 
        if(tabulist.size() > max_tabu_size) tabulist.pop_front() ; 
    }
    p = best ;
    return ;
}


/* Cutting Procedures */
/* Kaibel Box 1 Cuts */
void box1cuts(lpProblem& prob, int beta)
{
    random_device rd ; 
    mt19937 gen(rd()) ; 
    uniform_int_distribution<> distrib(1, n) ;
    uniform_int_distribution<> qdis(3, n-5) ; 

    GRBVar* vars = prob.model.getVars() ; 
    int cuts = 0, cnt_limit = 3000, i=0 ;
    double limit=1.0 ;  

    if(beta==3) limit=3.0 ; 

    omp_set_num_threads(omp_get_max_threads()) ; 
    #pragma omp parallel
    {
        while(i<cnt_limit)
        {   
            i++ ; 
            vector<int> p, q ; 
            vector<int> v(n) ;
            for(int i=0;i<n;i++) v[i] = i+1 ; 

            /* Initial random q  */
            shuffle(v.begin(), v.end(), rd) ; 
            int qsize = qdis(gen) ; 
            for(int i=0;i<qsize;i++) q.push_back(v[i]) ;

            /* Initial random p */
            shuffle(v.begin(), v.end(), rd) ; 
            int psize = qdis(gen) ;  /* For now, make it a constant */
            for(int i=0;i<psize;i++) p.push_back(v[i]) ; 
 
            if(psize + qsize > n-3+beta){i--; continue ; }
            /* cout << "P " ; 
            for(int i : p) cout << i << " " ; 
            cout << endl ;  */ 
            /* tabu search for approximately optimal solution */
            if(beta==2 && fitness(prob, p, q, 0, beta) < 0.6) { i--; continue ;}
            if(beta==3 && fitness(prob, p, q, 0, beta) < 2.5) {i--; continue; }
            //cout << "Before " << fitness(prob, p, q, 0, beta) << endl  ; 
            tabu_search(prob, p,q, psize, 0, beta) ; 
            //cout << "After " << fitness(prob, p, q, 0, beta) << endl << endl ; 
            /* Check limit  */
            if(fitness(prob, p, q, 0, beta) > limit)
            {
                cuts++ ; 
                GRBLinExpr exp ; 
                double pos= (double) beta-1, neg=-1.0 ; 
                for(int i : p)
                {
                    for(int j : q)
                    {
                        exp.addTerms(&pos, &vars[oneIndex(i,j,i,j)], 1) ; 

                        for(int k : p)
                        {
                            if(i>=k) continue ; 
                            for(int l : q) exp.addTerms(&neg, &vars[oneIndex(i,j,k,l)], 1) ; 
                        }   
                    }
                }

                prob.model.addConstr(exp, GRB_LESS_EQUAL, limit) ; 
            } 
        }
    }
    cout << cuts << endl ; 
    delete vars ; 
}

/* QAP2 Cuts */
void qap2(lpProblem& prob)
{
    int cuts = 0 ; 
    GRBVar* vars = prob.model.getVars() ; 

    for(int k=1;k<=n;k++)
    {
        for(int l=1;l<=n;l++)
        {
            for(int m=3;m<=3;m++)
            {
                vector<pair<double, pair<int, int>>> lst;

                for(int i=1;i<=n;i++)
                    for(int j=1;j<=n;j++) lst.push_back({prob.vars[oneIndex(i,j,k,l)], {i,j}}) ; 
                
                sort(lst.begin(), lst.end(), greater<>()) ;
            

            }
        }
    }

    cout << "Cuts Added : " << cuts << endl;
    delete vars ;  
    return  ; 
}

/* Branch and bound method */
void branchAndBound(lpProblem& prob)
{
    queue<lpProblem> q  ; 
    q.push(prob) ; 
    double expected_gap ; 

    int it = 0 ;  
    while(!q.empty())
    {
        if(it >5 ) break ; 
        it++ ; 
        cout << it << " Global Upper Bound : " << globalUpperBound << endl ;  

        if(lpRelaxation(q.front())) {q.pop(); continue;} 
        findUpperBound(q.front()) ; 

        cout << "Before adding cuts : " << q.front().lowerBound << " " << q.front().upperBound << endl ; 

        globalUpperBound = min(globalUpperBound, q.front().upperBound)  ;
        
        // Checking before adding cuts  
        if(q.front().lowerBound >= globalUpperBound || double_equal(q.front().lowerBound, q.front().upperBound)){ q.pop() ; continue ; }

        expected_gap = q.front().lowerBound/globalUpperBound ; 
        //  Cutting 
        /* for(int i=1;i<=1;i++){
            if(i & 1 && expected_gap < 0.975) box1cuts(q.front()) ; 
            else qap2(q.front()) ;  

            // Re-solving after adding Cuts 
            lpRelaxation(q.front()) ; 
            cout << q.front().lowerBound << endl ; 
            findUpperBound(q.front()) ; 

            globalUpperBound = min(globalUpperBound, q.front().upperBound)  ;
        } */

        // Checking again after adding cuts  
        if(q.front().lowerBound >= globalUpperBound || double_equal(q.front().lowerBound, q.front().upperBound)){ q.pop() ; continue ; }

        cout << "After adding cuts : " << q.front().lowerBound << " " << q.front().upperBound << endl ; 
         

        // Branching on max value found on diagonal 
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
        if(p.first==0 || p.second==0){q.pop(); continue;}

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
    //branchAndBound(prob) ;
    lpRelaxation(prob) ; 
    double last = 0.0, current = prob.lowerBound; 
    int i = 1 ; 
    cout <<"Initial Solution : " << current << endl ;
    cout << "Number of threads availabe " << omp_get_max_threads() << endl ; 
    while(abs(last - current)>0.1)
    {
        auto t_start = chrono::steady_clock::now() ; 
        cout << "Iteration : " << i << endl ; 
        qap2(prob) ; 
        last = current ;
        auto t_end = chrono::steady_clock::now() ; 
        lpRelaxation(prob) ; 
        current = prob.lowerBound ; 
        cout << "Solution " << current << endl ;
        cout << "Time Taken to find cuts :  " << (double) chrono::duration<double, milli>(t_end - t_start).count()/60000 << "min" << endl << endl ; 
        i++ ; 
    }
    last = 0.0, i=1 ; 
    /* while(abs(last -current)>0.1)
    {
        auto t_start = chrono::steady_clock::now() ; 
        cout << "Iteration : " << i << endl ; 
        box1cuts(prob, 2) ;
        last = current ;  
        auto t_end = chrono::steady_clock::now() ; 
        lpRelaxation(prob) ; 
        current = prob.lowerBound ; 
        cout << "Solution " << current << endl ;
        cout << "Time Taken to find cuts :  " << (double) chrono::duration<double, milli>(t_end - t_start).count()/60000 << "min" << endl << endl ; 
        i++ ;  
    } 
    last = 0.0, i=1;  
    while(abs(last-current)>0.1)
    {
        auto t_start = chrono::steady_clock::now() ; 
        cout << "Iteration : " << i << endl ; 
        box1cuts(prob, 3) ; 
        last = current ; 
        auto t_end = chrono::steady_clock::now() ; 
        lpRelaxation(prob) ; 
        current = prob.lowerBound ; 
        cout << "Solution " << current << endl;
        cout << "Time Taken to find cuts :  " << (double) chrono::duration<double, milli>(t_end - t_start).count()/60000 << "min" << endl << endl ; 
        i++ ; 
    }  
 */
    /* for(int i=0;i<10;i++)
    {
        cout << "Iteration " << i << endl ; 
        if(!i) lpRelaxation(prob) ; 
        cout << "Solution " << prob.lowerBound << endl ; 
        box1cuts(prob);
        lpRelaxation(prob) ;   
        cout << "Solution after adding cuts " << prob.lowerBound << endl << endl << endl ;
    } */ 

    auto end = chrono::steady_clock::now() ; 
    cout << "Time Taken :  " << (double) chrono::duration<double, milli>(end - start).count()/60000 << "min" << endl ; 

    return 0 ; 
}