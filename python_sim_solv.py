from scipy import optimize
import numpy as np

print("------------------------------------------------------------")

ob_func = np.array([-1,-9,-1])
ub_coeff_matrix = np.array([[1,2,3],
                   [3,2,2]
                   ])
ub_rhs_matrix = np.array([9,15])

bounds = (0, None) # already default params, but specifying just for show
soln = optimize.linprog(c=ob_func,A_ub=ub_coeff_matrix,b_ub=ub_rhs_matrix, method='simplex', bounds=bounds)
print("Project Solve")
print(soln)

# Should assert:

'''
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

'''
