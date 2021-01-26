#include<stdio.h>
#include<stdlib.h>
#include<glpk.h>
#include<string.h>


int main(){
    glp_prob *lp ; 
    glp_tran *tran ; 
    int ret ; 

    lp = glp_create_prob() ; 
    tran = glp_mpl_alloc_wksp() ; 
    ret = glp_mpl_read_model(tran, "model.mod", 0) ; 
    if(ret){
        fprintf(stderr, "Error on translating model.\n") ; 
        goto skip ;
    }
    /*ret = glp_mpl_read_data(tran, "data.dat") ; 
    if(ret){
        fprintf(stderr, "Error on reading data.\n") ; 
        goto skip ; 
    }*/
    ret = glp_mpl_generate(tran, NULL) ; 
    if(ret){
        fprintf(stderr, "Error on generating mode.\n") ; 
        goto skip ; 
    }
    glp_mpl_build_prob(tran, lp) ; 
    glp_create_index(lp) ;
    
    char *fname = "model.mps" ;
    //ret = glp_write_lp(lp,NULL,fname) ;  
    ret = glp_write_mps(lp, GLP_MPS_FILE, NULL, fname) ; 
    if(ret){
        fprintf(stderr, "Error while generating mps file.\n") ; 
    }
    /*
    glp_smcp param ; 
    glp_init_smcp(&param) ;
    param.presolve = GLP_ON ;  

    glp_create_index(lp) ; 
    ret = glp_simplex(lp, &param) ; 

    glp_print_sol(lp, "sol.txt") ;  
    */

    skip : 
        glp_mpl_free_wksp(tran) ; 
        glp_delete_prob(lp) ; 
        return 0 ; 
}