#include<bits/stdc++.h>
#include<glpk.h> 
#include<filesystem>

using namespace std ;
using recursive_directory_iterator = std::filesystem::recursive_directory_iterator ; 

int main()
{
    string path ="qapins" ; 

    for(const auto& file : recursive_directory_iterator(path))
    {
        cout << file << endl ;  
        ofstream outfile ; 
        outfile.open("model.mod", ios_base::out) ; 
        ifstream infile_tmp, infile ; 
        string news = file.path() ;
        string output_name = news.substr(7, news.length()-10) ; 
        output_name = output_name.append("mps") ;  
        infile.open(news, ios_base::in) ;

        infile_tmp.open("tmp", ios_base::in) ; 
        outfile << infile_tmp.rdbuf() ;  

        outfile << endl << "param n := " ; 
        
        int n ; 
        infile >> n ; 
        outfile << n << " ;" << endl << endl << "param a : " ; 

        for(int i=1;i<=n;i++)   
            outfile << i << " " ;
        outfile << ":=" << endl ; 

        for(int i=1;i<=n;i++)
        {
            outfile << i << " " ; 
            int tmp ; 
            for(int j=1;j<=n;j++)
            {
                infile >> tmp ; 
                outfile << tmp << " " ; 
            }
            if(i!=n) outfile << endl ; 
            else outfile << ";" << endl ;   
        }

        outfile << endl << endl << "param b : " ; 

        for(int i=1;i<=n;i++)   
            outfile << i << " " ;
        outfile << ":=" << endl ;

        for(int i=1;i<=n;i++)
        {
            outfile << i << " " ; 
            int tmp ; 
            for(int j=1;j<=n;j++)
            {
                infile >> tmp ; 
                outfile << tmp << " " ; 
            }
            if(i!=n) outfile << endl ; 
            else outfile << ";" << endl ;   
        }

        glp_prob *lp ; 
        glp_tran *tran ; 
        int ret ; 

        lp = glp_create_prob() ; 
        tran = glp_mpl_alloc_wksp() ; 
        ret = glp_mpl_read_model(tran, "model.mod", 0) ; 
        if(ret){
            fprintf(stderr, "Error on translating model.\n") ; 
            break ; 
        }
        ret = glp_mpl_generate(tran, NULL) ; 
        if(ret){
            fprintf(stderr, "Error on generating mode.\n") ; 
            break ;  
        }
        glp_mpl_build_prob(tran, lp) ; 
        glp_create_index(lp) ;
    
        char fname[output_name.length()] ;
        strcpy(fname, output_name.c_str()) ; 
        ret = glp_write_mps(lp, GLP_MPS_FILE, NULL, fname) ; 
        if(ret){
            fprintf(stderr, "Error while generating mps file.\n") ; 
        }
     
        glp_mpl_free_wksp(tran) ; 
        glp_delete_prob(lp) ; 
    }
    
    return 0 ;
}