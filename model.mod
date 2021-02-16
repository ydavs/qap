
param n, integer ; 

set I := 1..n ; 
set J := 1..n ; 
set K := 1..n ; 
set L := 1..n ; 

param a{i in I, j in J}, integer ; 
param b{i in I, j in J}, integer ; 

var y{i in I, j in J, k in K, l in L}, >=0; 

minimize cost: sum{i in I, j in J, k in K, l in L} a[i,k]*b[j,l]*y[i,j,k,l] ; 
s.t. ra{i in I, j in J, k in K, l in L}: y[i,j,k,l] = y[k,l,i,j] ; 
s.t. rb{i in I, j in J, l in L:j!=l}: y[i,j,i,l] = 0 ; 
s.t. rc{i in I, j in J, l in L:j!=l}: y[j,i,l,i] = 0 ; 
s.t. rd{i in I, j in J, l in L}: sum{k in K} y[i, j, k, l] = y[i,j,i,j] ; 
s.t. re{i in I, j in J, l in L}: sum{k in K} y[i, j, l, k] = y[i,j,i,j] ;
s.t. rf{i in I}: sum{j in J} y[i,j,i,j] = 1 ;
s.t. rg{i in I}: sum{j in J} y[j,i,j,i] = 1 ; 

data ; 

param n := 14;

param a: 1 2 3 4 5 6 7 8 9 10 11 12 13 14:= 
1  0  1  2  2  3  4  4  5  3  5  6  7  8  9
2  1  0  1  1  2  3  3  4  2  4  5  6  7  8
3  2  1  0  2  1  2  2  3  1  3  4  5  6  7
4  2  1  2  0  1  2  2  3  3  3  4  5  6  7
5  3  2  1  1  0  1  1  2  2  2  3  4  5  6
6  4  3  2  2  1  0  2  3  3  1  2  3  4  5
7  4  3  2  2  1  2  0  1  3  1  2  3  4  5
8  5  4  3  3  2  3  1  0  4  2  1  2  3  4
9  3  2  1  3  2  3  3  4  0  4  5  6  7  8
10  5  4  3  3  2  1  1  2  4  0  1  2  3  4
11  6  5  4  4  3  2  2  1  5  1  0  1  2  3
12  7  6  5  5  4  3  3  2  6  2  1  0  1  2
13  8  7  6  6  5  4  4  3  7  3  2  1  0  1
14  9  8  7  7  6  5  5  4  8  4  3  2  1  0;
 
param b: 1 2 3 4 5 6 7 8 9 10 11 12 13 14:=
1  0  3  4  6  8  5  6  6  5  1  4  6  1  5
2  3  0  6  3  7  9  9  2  2  7  4  7  9  6
3  4  6  0  2  6  4  4  4  2  6  3  6  5  6
4  6  3  2  0  5  5  3  3  9  4  3  6  3  4
5  8  7  6  5  0  4  3  4  5  7  6  7  7  3
6  5  9  4  5  4  0  8  5  5  5  7  5  1  8
7  6  9  4  3  3  8  0  6  8  4  6  7  1  8
8  6  2  4  3  4  5  6  0  1  5  5  3  7  5
9  5  2  2  9  5  5  8  1  0  4  5  2  4  5
10  1  7  6  4  7  5  4  5  4  0  7  7  5  6
11  4  4  3  3  6  7  6  5  5  7  0  9  6  5
12  6  7  6  6  7  5  7  3  2  7  9  0  6  5
13  1  9  5  3  7  1  1  7  4  5  6  6  0  5
14  5  6  6  4  3  8  8  5  5  6  5  5  5  0;