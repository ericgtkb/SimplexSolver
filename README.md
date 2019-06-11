# Simplex Solver
An implementation of the simplex algorithm for solving linear programs

## Input format
A Linear program in satndard form.
The first line specifies the number of constraints n and the number of variables m.
Each of the next n lines specifies the coefficients of the linear inequalities.
The next line contains n numbers for the right hand sides of the inequalities.
The last line specifies the coefficients of the objective function.

For example, the input for the linear program  
Maximize  
-x + 2y  
Subject to  
-x -y <= -1  
x <= 2  
y <= 2  
is  
```
3 2
-1 -1
1 0
0 1
-1 2 2
-1 2
```
