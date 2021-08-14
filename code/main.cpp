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
        prob.model.set(GRB_DoubleParam_IterationLimit, 50000) ; 
        //prob.model.set(GRB_IntParam_ScaleFlag, 2) ; 
        //prob.model.set(GRB_DoubleParam_ObjScale, -0.5) ; 
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
        else{ cout << "Optmization is stopped with status" << optimstatus << endl ; return -1; }
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

    for(int i=1;i<=n;i++)
    {
        for(int j=1;j<=n;j++) val+= a[i-1][j-1]*b[myJ[i]-1][myJ[j]-1] ; 
    }

    prob.upperBound = val ; 
    return prob.upperBound ; 
}

/* tabu Search */
double fitness(lpProblem &prob, int *p, int *q, int beta, int psize, int qsize)
{
     
    GRBVar* vars = prob.model.getVars() ; 
    double val=0.0 ; 

    for(int i=0;i<psize;i++)
    {
        for(int j=0;j<qsize;j++)
        {
            val+= (beta-1)*prob.vars[oneIndex(p[i], q[j], p[i], q[j])] ; 

            for(int k=0;k<psize;k++)
            {
                if(p[i]>=p[k]) continue ; 
                for(int l=0;l<qsize;l++)
                {
                    val -= prob.vars[oneIndex(p[i], q[j], p[k], q[l])] ; 
                }

            }
        }
    }

    delete vars ; 
    return val ;
}

bool chk_arr(int *p, int *q, int psize)
{
    int chk = 1 ; 
    for(int i=0;i<psize;i++)
    {
        if(p[i] != q[i]) chk =0 ; 
    }
    return chk ; 
}

void tabu_search(lpProblem &prob, int *p, int *q,  int psize, int qsize, int beta)
{  
    int best_candidate[psize] ;
    
    for(int i=0;i<psize;i++) best_candidate[i] = p[i] ;  
    
    int tabulist[max_tabu_size][psize] ; 
    
    for(int i=0;i<psize;i++) tabulist[0][i] = p[i] ; 

    int it = 0 ; 

    while(it < num_iteration) 
    {
        it++ ; 
        
        int valid[n-psize], bitarray[n+1], cbc[psize] ;
        
        for(int i=0;i<=n;i++) bitarray[i] = 0 ; 
        
        for(int i=0;i<psize;i++) 
        {
            bitarray[best_candidate[i]] = 1 ;
            cbc[i] = best_candidate[i] ; 
        }  
        
        for(int i=1,j=0;i<=n;i++)
        {
            if(bitarray[i] != 1)
            {
                valid[j] = i ; 
                j++ ; 
            }
        }

        int ft = 1 ; 

        for(int i=0;i<psize;i++)
        {
            int tmp[psize] ; 

            for(int j=0;j<psize;j++)
            {
                if(j==i) continue ; 
                tmp[j] = cbc[j] ; 
            }

            for(int j=0;j<n-psize;j++)
            {
                tmp[i] = valid[j] ; 

                if(ft)
                {
                    for(int k=0;k<psize;k++) best_candidate[k] = tmp[k] ; 
                    ft = 0 ; 
                }

                int is_eq = 0 ; 
                for(int k=0;k<max_tabu_size;k++) if(chk_arr(tmp, tabulist[k], psize)) is_eq = 1 ; 

                if(is_eq) continue ; 

                if(fitness(prob, tmp, q, beta, psize, qsize) > fitness(prob, best_candidate, q,  beta, psize, qsize))
                {
                    for(int k=0;k<psize;k++) best_candidate[k] = tmp[k] ; 

                    if(fitness(prob, best_candidate, q, beta, psize, qsize) > fitness(prob, p, q, beta, psize, qsize))
                    {
                        for(int k=0;k<psize;k++) p[k] = best_candidate[k] ;
                    }
                }
            }
        }
    }

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
    int cuts = 0, cnt_limit = 8000, i=0 ;
    double limit=1.0 ;  

    if(beta==3) limit=3.0 ; 
    
    omp_set_num_threads(omp_get_max_threads()) ; 
    #pragma omp parallel for shared(prob)  
    for(int i =0;i<cnt_limit;i++)
    {
        int v[n] ; 
        for(int i=0;i<n;i++) v[i] = i+1 ; 
          
        /* Initial random q  */
        random_shuffle(v, v + n) ; 
        int qsize = qdis(gen) ; 
        int q[qsize] ; 
        for(int i=0;i<qsize;i++) q[i] = v[i] ; 

        /* Initial random p */
        random_shuffle(v, v+n) ;  
        int psize = qdis(gen) ; 
        int p[psize] ; 
        for(int i=0;i<psize;i++) p[i] = v[i] ;  
 
        if(psize + qsize > n-3+beta){i--; continue ; }
        
        /* tabu search to improve cut */
        if(beta==2 && fitness(prob, p, q,beta, psize, qsize) < 0.6) { i--; continue ;}
        if(beta==3 && fitness(prob, p, q, beta, psize, qsize) < 2.5) {i--; continue; } 
        
        tabu_search(prob, p,q, psize, qsize, beta) ; 
        
        /* Add cut to the model */
        if(fitness(prob, p, q,beta, psize, qsize) > limit)
        {
            cuts++ ; 
            GRBLinExpr exp ; 
            double pos= (double) beta-1, neg=-1.0 ; 
            for(int i=0;i<psize;i++)
            {
                for(int j=0;j<qsize;j++)
                {
                    exp.addTerms(&pos, &vars[oneIndex(p[i],q[j],p[i],q[j])], 1) ; 

                    for(int k=0;k<psize;k++)
                    {
                        if(p[i]>=p[k]) continue ; 
                        for(int l=0;l<qsize;l++) exp.addTerms(&neg, &vars[oneIndex(p[i],q[j],p[k],q[l])], 1) ; 
                    }   
                }
            }

            #pragma omp critical
            prob.model.addConstr(exp, GRB_LESS_EQUAL, limit) ; 
        } 
    }
    cout << cuts << endl ; 
    delete vars ; 
}

/* QAP1 Cuts */
void qap1(lpProblem& prob)
{
    int cuts =0 ; 
    GRBVar* vars = prob.model.getVars() ; 

    for(int k=1;k<=n;k++)
    {
        for(int l=1;l<=n;l++)
        {
            vector<pair<double, pair<int,int>>> lst ; 

            for(int i=1;i<=n;i++)
                for(int j=1;j<=n;j++) lst.push_back({prob.vars[oneIndex(i,j,k,l)], {i,j}}) ;

            sort(lst.begin(), lst.end(), greater<>()) ; 

            for(int a=0;a<n;a++)
            {
                if(lst[a].second.first ==k || lst[a].second.second ==l ) continue ; 

                for(int b=0;b<n;b++)
                {
                    if(lst[b].second.first ==k || lst[b].second.second ==l || lst[a].second.first == lst[b].second.first || lst[a].second.second == lst[b].second.second) continue ; 

                    GRBLinExpr exp ; 
                    double pos=1.0, neg=-1.0, sum=0.0 ;
                    int p1=lst[a].second.first, q1=lst[a].second.second, p2=lst[b].second.first, q2=lst[b].second.second ;
                    exp.addTerms(&pos, &vars[oneIndex(p1,q1,k,l)],1) ; 
                    sum += prob.vars[oneIndex(p1,q1,k,l)] ; 
                    exp.addTerms(&pos, &vars[oneIndex(p2,q2,k,l)],1) ; 
                    sum += prob.vars[oneIndex(p2,q2,k,l)] ; 
                    exp.addTerms(&pos, &vars[oneIndex(p1,q2,k,l)],1) ; 
                    sum += prob.vars[oneIndex(p1,q2,k,l)] ; 
                    exp.addTerms(&neg, &vars[oneIndex(k,l,k,l)],1) ; 
                    sum -= prob.vars[oneIndex(k,l,k,l)] ; 
                    exp.addTerms(&neg, &vars[oneIndex(p1,q1,p2,q2)],1) ; 
                    sum -= prob.vars[oneIndex(p1,q1,p2,q2)] ; 

                    if(sum > 0.05)
                    {
                        prob.model.addConstr(exp, GRB_LESS_EQUAL, 0) ;
                        cuts++ ;  
                    }
                }
            }
        }
    }

    cout <<"Cuts Added : " << cuts << endl ; 
    return ; 
}

/* QAP2 Cuts */
void qap2(lpProblem& prob)
{
    int m = 3 ;
    int cuts = 0 ; 
    GRBVar* vars = prob.model.getVars() ; 

    for(int k=1;k<=n;k++)
    {
        for(int l=1;l<=n;l++)
        {
            vector<pair<double, pair<int, int>>> lst ; 
            
            for(int i=1;i<=n;i++)
                for(int j=1;j<=n;j++) lst.push_back({prob.vars[oneIndex(i,j,k,l)], {i,j}}) ; 

            sort(lst.begin(), lst.end(), greater<>()) ;

            for(int a=0;a<n;a++)
            {
                if(lst[a].second.first ==k || lst[a].second.second ==l ) continue ; 
                for(int b=0;b<n;b++)
                {
                    if(lst[b].second.first ==k || lst[b].second.second ==l || lst[a].second.first == lst[b].second.first || lst[a].second.second == lst[b].second.second) continue ; 
                    for(int c=0;c<n;c++)
                    {
                        if(lst[c].second.first ==k || lst[c].second.second ==l || lst[c].second.first == lst[b].second.first || lst[c].second.second == lst[b].second.second || lst[a].second.first == lst[c].second.first || lst[a].second.second == lst[c].second.second) continue ; 

                        vector<int> ivec, jvec ; 
                        ivec.push_back(lst[a].second.first) ; 
                        ivec.push_back(lst[b].second.first) ;
                        ivec.push_back(lst[c].second.first) ;

                        jvec.push_back(lst[a].second.second) ; 
                        jvec.push_back(lst[b].second.second) ;
                        jvec.push_back(lst[c].second.second) ;

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
                            cuts++ ;  
                        }
                    }

                    vector<int> ivec, jvec ; 
                    ivec.push_back(lst[a].second.first) ; 
                    ivec.push_back(lst[b].second.first) ; 

                    jvec.push_back(lst[a].second.second) ; 
                    jvec.push_back(lst[b].second.second) ;

                    GRBLinExpr exp ; 
                    double val = 0.0, tmp1=1.0, tmp2=-1.0 ; 
                    for(int i=0;i<2;i++)
                    {
                        val += prob.vars[oneIndex(ivec[i], jvec[i], k, l)] ; 
                        exp.addTerms(&tmp1, &vars[oneIndex(ivec[i], jvec[i], k, l)], 1) ; 
                    } 

                    val+= prob.vars[oneIndex(ivec[0], jvec[1], k, l)] - prob.vars[oneIndex(k, l, k, l)] - prob.vars[oneIndex(ivec[1], jvec[1], ivec[0], jvec[0])] ; 
                    exp.addTerms(&tmp1, &vars[oneIndex(ivec[0], jvec[1], k, l)], 1) ; 
                    exp.addTerms(&tmp2, &vars[oneIndex(k, l, k, l)], 1) ; 
                    exp.addTerms(&tmp2, &vars[oneIndex(ivec[1], jvec[1], ivec[0], jvec[0])], 1) ; 
                    
                    if(val > 0.05)
                    {
                        prob.model.addConstr(exp, GRB_LESS_EQUAL, 0) ; 
                        cuts++ ; 
                    }
                }
            }
        }
    } 

    cout << "Cuts Added : " << cuts << endl; 
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
        if(it >100) break ; 
        it++ ; 
        cout << it << " Global Upper Bound : " << globalUpperBound << endl ;  

        if(lpRelaxation(q.front())) {q.pop(); continue;} 
        findUpperBound(q.front()) ; 

        cout << "Stats : " << q.front().lowerBound << " " << q.front().upperBound << endl ; 

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
    cout.precision(5) ; 

    if(argc < 3)
    {
        cout << "Usage : ./main model instance" << endl ; 
        return 1;  
    }

    readCoeff(argv[2]) ; 
    lpProblem prob = initialize(argv[1]) ; 

    int is_bbc = 0 ; 

    if(is_bbc) 
    {
        /* Branch and Cut */   
        branchAndBound(prob) ;
    }
    else
    {
        /* Plane Cutting */
        lpRelaxation(prob) ; 
        cout <<"Initial Solution : " << prob.lowerBound << endl ;
        cout << "Number of threads availabe " << omp_get_max_threads() << endl ; 

        double lambda = 2 ; 
        for(int j=3;j<=4;j++)
        {
            double last = 0.0, current = prob.lowerBound ; 
            int itr = 1 ; 
            while(abs(last - current) > lambda && itr<14)
            {
		if(j<3 && itr>2) break ; 
                auto tcc_start = chrono::steady_clock::now() ; 
                cout << "Cut " << j << " Iteration : " << itr << endl ;
                switch(j)
                {
                    case 1: 
                        qap1(prob) ;
                        break ; 
                    case 2:
                        qap2(prob) ; 
                        break ; 
                    case 3:
                        box1cuts(prob, 2) ; 
                        break ; 
                    case 4:
                        box1cuts(prob, 3) ; 
                        break ; 
                } 
                auto tcc_end =  chrono::steady_clock::now() ; 
                auto tlp_start = chrono::steady_clock::now() ; 
                int x = lpRelaxation(prob) ;
                if(x==-1) continue ;  
                last = current ;
                current =  prob.lowerBound ; 
                auto tlp_end = chrono::steady_clock::now() ; 
                cout << "Solution " << current << "Time for cuts : " << (double) chrono::duration<double, milli>(tcc_end - tcc_start).count()/60000 << "Min" << " Time for LP : " << (double) chrono::duration<double, milli>(tlp_end-tlp_start).count()/60000 <<"Min" << endl ; 
                itr++ ; 
            }

        }
        prob.model.write("lipa20a_bb.mps") ;
    }


    auto end = chrono::steady_clock::now() ; 
    cout << "Time Taken :  " << (double) chrono::duration<double, milli>(end - start).count()/60000 << "min" << endl ; 

    return 0 ; 
}
