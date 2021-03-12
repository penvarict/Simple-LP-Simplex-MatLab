# Project 1 	

**By: Charlie Penvari **

The output for my implemented solution in MATLAB 2019b is as follows: 

```
Project problem
------------------------------------------------
Basic variables: x_4 5  Non-basic variables: x_1 2 3 
Entering basic variable: x_2 Leaving basic variable: x_4
Current cpf 0.00

Basic variables: x_2 5  Non-basic variables: x_1 4 3 
Entering basic variable: x_4 Leaving basic variable: x_2
Current cpf 40.50

End simplex loop
Optimal value is Z = 40.500
The value occurs at the variable x_2 = 4.500
The value occurs at the variable x_5 = 6.000
All other variables are presumed to equal zero
```



We also do the Wyndor problem just to show that it works: 

```
Wyndor Problem
------------------------------------------------
Basic variables: x_3 4 5  Non-basic variables: x_1 2 
Entering basic variable: x_2 Leaving basic variable: x_4
Current cpf 0.00

Basic variables: x_3 2 5  Non-basic variables: x_1 4 
Entering basic variable: x_1 Leaving basic variable: x_5
Current cpf 30.00

Basic variables: x_3 2 1  Non-basic variables: x_5 4 
Entering basic variable: x_5 Leaving basic variable: x_1
Current cpf 36.00

End simplex loop
Optimal value is Z = 36.000
The value occurs at the variable x_3 = 2.000
The value occurs at the variable x_2 = 6.000
The value occurs at the variable x_1 = 2.000
All other variables are presumed to equal zero
```

In the solver from scipy.optimize we get: 

```
------------------------------------------------------------
Project Solve
     con: array([], dtype=float64)
     fun: -40.5
 message: 'Optimization terminated successfully.'
     nit: 4
   slack: array([0., 6.])
  status: 0
 success: True
       x: array([0. , 4.5, 0. ])

```

The above output corresponds with the MATLAB script. 