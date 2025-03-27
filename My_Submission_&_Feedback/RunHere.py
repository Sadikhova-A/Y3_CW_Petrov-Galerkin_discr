"""
MATH3036 Coursework 1 main script

@author: Kris van der Zee (Lecturer)

version: Jan 2025
"""



#%% Question 1 (Compute PG approximation for N=2)

import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as DE


# Set problem data
# Function in the differential equation
def f(x):
    return np.pi/2 * np.cos(np.pi/2*x)

# Number of subintervals to be used in numerical integration
n = 100

# Assemble the 2x2 system
A = DE.AssemblePGmatrixN2()
b = DE.AssemblePGvectorN2(f,n)

print("A =",A)
print("b =",b)

# Solve system
u = np.linalg.solve(A,b)
print("u =",u)


# Plot approximation
N_p = 100
x_p = 1/N_p * np.arange(0,N_p+1,1)
u_p = u[0]*x_p + u[1]/2*(x_p**2)   # Approximation at plot points
fig, ax = plt.subplots()  # Create a figure containing a single axes.
ax.plot(x_p,u_p,'b-')

# Exact solution
def u_ex(x):
    return np.sin(np.pi/2*x)
# Plot exact solution
u_ex_p = u_ex(x_p)    # Exact solution at plot points
ax.plot(x_p,u_ex_p,'r-')

# Legend
ax.legend(['u_h(x)','u(x)'])



#%% Question 2 (Compute PG approximation for N=3)

import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as DE


# Set problem data
# Function in the differential equation
def f(x):
    return np.pi/2 * np.cos(np.pi/2*x)

# Number of subintervals to be used in numerical integration
n = 100

# Assemble the 2x2 system
A = DE.AssemblePGmatrixN3()
b = DE.AssemblePGvectorN3(f,n)

print("A =",A)
print("b =",b)

# Solve system
u = np.linalg.solve(A,b)
print("u =",u)

# Plot approximation
N_p = 100
x_p = 1/N_p * np.arange(0,N_p+1,1)
u_p = u[0]*x_p + u[1]/2*(x_p**2) + u[2]/3*(x_p**3)   # Approximation at plot points
fig, ax = plt.subplots()  # Create a figure containing a single axes.
ax.plot(x_p,u_p,'b-')

# Exact solutuion
def u_ex(x):
    return np.sin(np.pi/2*x)
# Plot exact solution
u_ex_p = u_ex(x_p)    # Exact solution at plot points
ax.plot(x_p,u_ex_p,'r-')

# Legend
ax.legend(['u_h(x)','u(x)'])



#%% Question 3 (Compute PG approximation for N)

import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as DE


# Set problem data
# Function in the differential equation
def f(x):
    return np.pi/2 * np.cos(np.pi/2*x)

# Number of subintervals to be used in numerical integration
n = 100

# Degree
N = 4

# Assemble the 2x2 system
A = DE.AssemblePGmatrix(N)
b = DE.AssemblePGvector(N,f,n)
print("A =",A)
print("b =",b)

# Solve system
u = np.linalg.solve(A,b)
print("u =",u)


# Plot approximation
N_p = 100
x_p = 1/N_p * np.arange(0,N_p+1,1)
u_p = np.zeros(N_p+1)     # Initialize approximation at plot points
for j in np.arange(N)+1:
    u_p = u_p + u[j-1]/j*(x_p**j)
fig, ax = plt.subplots()  # Create a figure containing a single axes.
ax.plot(x_p,u_p,'b-')

# Exact solutuion
def u_ex(x):
    return np.sin(np.pi/2*x)
# Plot exact solution
u_ex_p = u_ex(x_p)    # Exact solution at plot points
ax.plot(x_p,u_ex_p,'r-')

# Legend
ax.legend(['u_h(x)','u(x)'])


#%% Question 4 (PG for General 1st-order differential equation)


import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as DE


# Set problem data
# Function in the differential equation
def f(x):
    return np.ones(np.shape(x))  # the function f(x) = 1

# Parameter in differential equation
c = 4

# Number of subintervals to be used in numerical integration
n = 100

# Degree
N = 3

# Assemble the 2x2 system
A = DE.AssemblePGmatrixGenDE(N,c)
b = DE.AssemblePGvector(N,f,n)
print("A =",A)
print("b =",b)

# Solve system
u = np.linalg.solve(A,b)
print("u =",u)


# Plot approximation
N_p = 100
x_p = 1/N_p * np.arange(0,N_p+1,1)
u_p = np.zeros(N_p+1)     # Initialize approximation at plot points
for j in np.arange(N)+1:
    u_p = u_p + u[j-1]/j*(x_p**j)
fig, ax = plt.subplots()  # Create a figure containing a single axes.
ax.plot(x_p,u_p,'b-')

# Exact solutuion
def u_ex(x):
    return -1/c * np.exp(-c*x) + 1/c
# Plot exact solution
u_ex_p = u_ex(x_p)    # Exact solution at plot points
ax.plot(x_p,u_ex_p,'r-')

# Legend
ax.legend(['u_h(x)','u(x)'])


#%% Question 5 (Galerkin FEM for 2nd-order problem)

import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as DE

# Test of vector b
# Set mesh parameter (N-1 is number of interior nodes)
N = 4;
# Set problem data
def f(x):
    return np.sin(np.pi*x)  # the function f(x) = 1
# Assemble system
b = DE.AssembleFEMvector(N,f)
print("b =",b)


# Test of matrix A and vector b
# Set mesh parameter (N-1 is number of interior nodes)
N = 4;

# Set problem data
def f(x):
    return np.ones(np.shape(x))  # the function f(x) = 1
c = 6

# Assemble system
A = DE.AssembleFEMmatrix(N,c)
print("A =",A)
b = DE.AssembleFEMvector(N,f)
print("b =",b)

# Solve system
u = np.linalg.solve(A,b)
print("u =",u)

# Plot FEM approximation
h = 1/N 
x_all = h*np.arange(0,N+1,1)
u_all = np.r_[0,u,0]  
fig, ax = plt.subplots()  # Create a figure containing a single axes.
ax.plot(x_all,u_all,'bo-')


# Exact solutuion
def u_ex(x):
    a = (np.cosh(np.sqrt(c)) - 1) /c * 1 / np.sinh(np.sqrt(c))
    b = -1/c
    return  a * np.sinh(np.sqrt(c)*x) + b * np.cosh(np.sqrt(c)*x) + 1/c

# Plot exact solution
N_p = 100
x_p = 1/N_p * np.arange(0,N_p+1,1)
u_ex_p = u_ex(x_p)
ax.plot(x_p,u_ex_p,'r-')

# Legend
ax.legend(['u_h(x)','u(x)'])






#%% Question 6 (Heat Equation using backward Euler - Galerkin FEM)


import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as DE

# Set mesh parameter 
alpha = 1/2
N = 8
tau = 0.02
Ntime = 16

# Set initial condition function
def u0_func(x):
    return 3*x*(1-x) + np.sin(2*np.pi*x)

# Approximation of u0
u0_vec = DE.AssembleVector_u0(N,u0_func)
print("u0_vec =",u0_vec)

# Matrices
M = DE.AssembleMatrix_M(N)
print("M =",M)
K = DE.AssembleMatrix_K(N)
print("K =",K)

# Solve heat equation
u_array = DE.HeatEqFEM(tau,alpha,Ntime,N,u0_func)
print("u_array =",u_array)


# Plot approximations
fig, ax = plt.subplots()  # Create a figure containing a single axes.
h = 1/N
x_all = h*np.arange(0,N+1,1)
for k in np.arange(Ntime+1):
    uk = u_array[k,:]
    uk_all = np.r_[0,uk,0]
    ax.plot(x_all,uk_all,'o-')





    

#%% Question 7  (Compute the 2d FEM system for N = 6)


import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as DE

# Assemble A when N=6
A = DE.Assemble2dFEMmatrixN6()
print("A = ",A)

# Assemble b
N = 6
b = DE.Assemble2dFEMvector(N)
print("b =",b)

# Solve system
u = np.linalg.solve(A,b)
print("u =",u)


#%% Question 8 (Compute the 2d FEM system for arbitrary N)


import numpy as np
import matplotlib.pyplot as plt
import DEtools_SUBMISSION as FEM

# Large N
N = 24
A = DE.Assemble2dFEMmatrix(N)
print("A = ",A)

b = DE.Assemble2dFEMvector(N)
print("b =",b)

# Solve system
u = np.linalg.solve(A,b)
print("u =",u)

# Plot the approximation using a scatter plot
fig = plt.figure()  
ax = fig.add_subplot(projection='3d')

# Plot boundary nodes
h = 1/N
xb = h*np.arange(0,N+1)
yb = 0*np.ones(N+1)
for j in np.arange(1,N+1):
    xb = np.hstack((xb,h*np.array([0,N-j])))
    yb = np.hstack((yb,j*h*np.array([1,1])))
zb = np.zeros(np.shape(xb))
ax.scatter(xb, yb, zb,c=zb,s=40)

# Plot interior nodes
h = 1/N
xs = []
ys = []
for j in np.arange(1,N):
    xs = np.hstack((xs,h*np.arange(1,N-j)))
    ys = np.hstack((ys,j*h*np.ones(N-1-j)))
zs = u
ax.scatter(xs, ys, zs,c=zs,s=40)

# Postprocess figure
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_zlim([min([0,min(u)]),max([0.03,max(u)])])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.view_init(elev=30., azim=-30)


