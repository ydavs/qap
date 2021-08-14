# qap
## Dependencies : 
* Gurobi (https://www.gurobi.com/documentation/9.1/remoteservices/linux_installation.html)
* Problem Instances(can be downloaded from http://anjos.mgi.polymtl.ca/qaplib/inst.html)

## Running Instructions : 

* Please make sure you've installed all the required dependencies
* Use make for compilation 

         make 
* Usage ./main {model file} {instance file} 

        ./main model.mps ins.dat 


Look at logfile.log for a detailed output of program. Output from QUBO SOLVER is not printed in logfile due to the large number of repetitions in a single run. 

When changing problem instances, make sure to update values in all, ins.dat, model.mod and main.h. 
