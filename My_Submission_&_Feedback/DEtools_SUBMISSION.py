"""
MATH3036 CW1 DEtools module

@author: Aylin Sadikhova
"""

import numpy as np



def AssemblePGmatrixN2():
    """
    Returns a numpy array of the 
    2x2 Petrov-Galerkin matrix for the bilinear form
       b(u,v) := \int_0^1 u' v dx 
    using trial space 
       Uh := {polynomials of degree at most 2, and being = 0 for x=0}, 
    with basis { x, x**2/2 } 
    and using the test space
       Vh := {polynomials of degree at most 1}
    with basis { 1, x }
                
    Parameters
    ----------
    none

    Returns
    -------
    A : numpy.ndarray, shape (2,2)
        Array containing the Petrov-Galerkin matrix.
    """
    
    # Assemble matrix
    A = np.array([[1,1/2],[1/2,1/3]])
        
    return A


def AssemblePGvectorN2(f,n):
    """
    Returns a numpy array of the 
    vector for the linear form
       l(v) := \int_0^1 f v dx 
    using test space
       Vh := {polynomials of degree at most 1}
    with basis { 1, x }.
    The integral is approximated numerically using
    a composite mid-point rule on n subintervals.
    
    Parameters
    ----------
    f : function
        Input function.
    n : integer
        Number of subintervals used in the 
        composite mid-point rule
        
    Returns
    -------
    b : numpy.ndarray, shape (2,)
        Array containing the vector.
    """
    
    # Define mid points (center of each subinterval)
    h = 1/n;
    xc = h*np.arange(n) + h/2

    # Integral approximated using composite mid-point rule 
    b1 = np.sum( h*f(xc) )
    b2 = np.sum( h*(f(xc)*xc) )
    
    # Assemble vector
    b = np.array([b1,b2])
    
    return b



def AssemblePGmatrixN3():
    """
    Returns a numpy array of the 
    3x3 Petrov-Galerkin matrix for the bilinear form
       b(u,v) := \int_0^1 u' v dx 
    using trial space 
       Uh := {polynomials of degree at most 3, and being = 0 for x=0}, 
    with basis { x, (x**2)/2, (x**3)/3 } 
    and using the test space
       Vh := {polynomials of degree at most 2}
    with basis { 1, x, x**2 }
                
    Parameters
    ----------
    none

    Returns
    -------
    A : numpy.ndarray, shape (3,3)
        Array containing the Petrov-Galerkin matrix.
    """
    
    # Assemble matrix
    A = [[1, 1/2, 1/3], [1/2, 1/3, 1/4], [1/3, 1/4, 1/5]]
        
    return A


def AssemblePGvectorN3(f,n):
    """
    Returns a numpy array of the 
    vector for the linear form
       l(v) := \int_0^1 f v dx 
    using test space
       Vh := {polynomials of degree at most 2}
    with basis { 1, x, x**2 }.
    The integral is approximated numerically using
    a composite mid-point rule on n subintervals.
    
    Parameters
    ----------
    f : function
        Input function.
    n : integer
        Number of subintervals used in the 
        composite mid-point rule
        
    Returns
    -------
    b : numpy.ndarray, shape (3,)
        Array containing the vector.
    """

    # Defining mid-points
    h = 1/n
    xc = h * np.arange(n) + h/2

    # Approximating integrals using the composite mid-point rule
    b1 = np.sum(h*f(xc))
    b2 = np.sum(h*f(xc)*xc)
    b3 = np.sum(h*f(xc)*xc**2)

    # Assembling the vector
    b = [b1, b2, b3]

    return b


def AssemblePGmatrix(N):
    """
    Returns a numpy array of the 
    NxN Petrov-Galerkin matrix for the bilinear form
       b(u,v) := \int_0^1 u' v dx 
    using trial space 
       Uh := {polynomials of degree at most N, and being = 0 for x=0}, 
    with basis { x, (x**2)/2, (x**3)/3,..., (x**N)/N } 
    and using the test space
       Vh := {polynomials of degree at most N-1}
    with basis { 1, x, x**2,..., x**(N-1) }
                
    Parameters
    ----------
    N : integer
        Degree of the polynomial approximation.

    Returns
    -------
    A : numpy.ndarray, shape (N,N)
        Array containing the Petrov-Galerkin matrix.
    """
    # Initiating a matrix
    A = np.zeros((N, N))

    # Computing the matrix
    for i in range(N):
        for j in range(N):
            A[i][j] = 1 / (i + j + 1)
            
    return A


def AssemblePGvector(N,f,n):
    """
    Returns a numpy array of the 
    vector for the linear form
       l(v) := \int_0^1 f v dx 
    using test space
       Vh := {polynomials of degree at most N-1}
    with basis { 1, x, x**2,..., x**(N-1) }.
    The integral is approximated numerically using
    a composite mid-point rule on n subintervals.
    
    Parameters
    ----------
    N : integer
        N-1 is the maximum polynomial degree of the test space.
    f : function
        Input function.
    n : integer
        Number of subintervals used in the 
        composite mid-point rule
        
    Returns
    -------
    b : numpy.ndarray, shape (N,)
        Array containing the vector.
    """
    # Defining mid-points
    h = 1/n
    xc = h * np.arange(n) + h/2

    # Initiating a vector
    b = np.zeros(N)

    # Assembling the vector
    for i in range(N):
        b[i] = np.sum(h * f(xc) * xc ** i)

    return b



def AssemblePGmatrixGenDE(N,c):
    """
    Returns a numpy array of the
    NxN Petrov-Galerkin matrix for the bilinear form
       b(u,v) := \int_0^1 u' v + c u v dx
    using trial space
       Uh := {polynomials of degree at most N, and being = 0 for x=0},
    with basis { x, (x**2)/2, (x**3)/3,..., (x**N)/N }
    and using the test space
       Vh := {polynomials of degree at most N-1}
    with basis { 1, x, x**2,..., x**(N-1) }

    Parameters
    ----------
    N : integer
        Degree of polynomial approximation.
    c : float
        Constant in front of the u v term in the PG integral.
    
    Returns
    -------
    A: numpy.ndarray, shape (N,N)
       Array containing the Petrov-Galerkin matrix.

    """
    # Initializing the matrix
    A1 = np.zeros((N, N))
    A2 = np.zeros((N, N))

    # Assemble matrix
    for i in range(N):
        for j in range(N):
            A1[i][j] = 1 / (i + j + 1)
            A2[i][j] = c / ((i + j + 2) * (j + 1))

    A = A1 + A2

    return A



def AssembleFEMmatrix(N,c):
    """
    Returns a numpy array of the
    (N-1)x(N-1) Galerkin FEM matrix for the bilinear form
       b(u,v) := \int_0^1 u' v' + c u v dx
    using trial space
       Uh := {continuous piecewise linears, equal to 0 at x = 0 and x = 1},
    with basis of the usual hat functions \phi_j(x_i) = 1 when i = j and 0 otherwise
    using FEM test space
       Vh := continuous piecewise linears,
             and being = 0 for x=0 and x=1,
    with basis { psi_i } being the hat function (node-wise),
    for a uniform mesh with N+1 grid-nodes.

    Parameters
    ----------
    N : integer
        N - 1 is the number of basis functions.
    c : float
        Non-negative constant in front of the u v term in the differential equation.
    
    Returns
    -------
    A: numpy.ndarray, shape (N-1,N-1)
        Array containing the Petrov-Galerkin matrix.

    """
    # Initializing width parameter
    h = 1/N

    # Initializing matrices whose linear combination equal the final matrix
    A1 = np.zeros((N-1, N-1))
    A2 = np.zeros((N - 1, N - 1))

    # Assembling the matrix
    for i in range(N-1):
        for j in range(N-1):
            if i == j:
                A1[i][j] = 2
                A2[i][j] = 4
            if i == j + 1 or j == i + 1:
                A1[i][j] = -1
                A2[i][j] = 1

    # Computing the matrix
    A = A1/h + c * h * A2 / 6

    return A





def AssembleFEMvector(N,f):
    """
    Returns a numpy array of the 
    FEM vector for the linear form
       l(v) := \int_0^1 f v dx 
    using FEM test space
       Vh := continuous piecewise linears,
             and being = 0 for x=0 and x=1, 
    with basis { psi_i } being the hat function (node-wise),
    for a uniform mesh with N+1 grid-nodes.  
    
    Parameters
    ----------
    N : integer
        N+1 is the number of nodes in the mesh.
    f : function
        Input function.
        
    Returns
    -------
    b : numpy.ndarray, shape (N-1,)
        Array containing the FEM vector.
    """
    # Initializing the mid-points
    h = 1 / N
    xm = h * (np.arange(1, N) - 1/2)
    xp = h * (np.arange(1, N) + 1/2)

    # Initializing the vector
    b = np.zeros(N-1)

    # Assembling the vector
    b = (f(xp) + f(xm)) * h / 2
    
    return b


def AssembleMatrix_M(N):

    # Initializing width parameter
    h = 1 / N

    # Initializing matrix
    M = np.zeros((N - 1, N - 1))  # (N-1)x(N-1) Dimensional Matrix

    # Assembling the matrix
    for i in range(N - 1):
        for j in range(N - 1):
            if i == j:
                M[i][j] = 4
            if i == j + 1 or j == i + 1:
                M[i][j] = 1

    # Scaling the matrix
    M = (h/6) * M

    return M

    
def AssembleMatrix_K(N):
    # Initializing width parameter
    h = 1 / N

    # Initializing matrices whose linear combination equal the final matrix
    K = np.zeros((N - 1, N - 1))  # (N-1)x(N-1) Dimensional Matrix

    # Assembling the matrix
    for i in range(N - 1):
        for j in range(N - 1):
            if i == j:
                K[i][j] = 2
            if i == j + 1 or j == i + 1:
                K[i][j] = -1

    # Computing the matrix
    K = K / h
    return K


def AssembleVector_u0(N, u0_func):
    # Initializing the mid-points
    h = 1 / N
    xm = h * (np.arange(1, N) - 1 / 2)
    xp = h * (np.arange(1, N) + 1 / 2)

    # Calling the M matrix
    M = AssembleMatrix_M(N)

    # Assembling the vector
    b = (u0_func(xp) + u0_func(xm)) * h / 2  # Same equation as in Q5

    # Computing u0
    u0_vec = np.linalg.solve(M, b)

    return u0_vec


def HeatEqFEM(tau, alpha, Ntime, N, u0_func):
    # Initializing the array
    u_array = np.zeros((Ntime+1, N-1))

    # Calling functions
    u_array[0] = AssembleVector_u0(N, u0_func)
    K = AssembleMatrix_K(N)
    M = AssembleMatrix_M(N)

    # Computing the u_array for i = 1, ..., Ntime + 1
    for i in range(1, Ntime + 1):
        c = np.matmul(M, u_array[i - 1])
        u_array[i] = np.linalg.solve((M + tau * alpha * K), c)
    
    return u_array




def Assemble2dFEMvector(N):
    # Setting step width
    h = 1/N

    # Initializing the vector
    dim = int((N - 2) * (N - 1) / 2)
    b = np.ones(dim) * h ** 2

    return b


def Assemble2dFEMmatrixN6():
    A = np.array([
        [4, -1, 0, 0, -1, 0, 0, 0, 0, 0],
        [-1, 4, -1, 0, 0, -1, 0, 0, 0, 0],
        [0, -1, 4, -1, 0, 0, -1, 0, 0, 0],
        [0, 0, -1, 4, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 4, -1, 0, -1, 0, 0],
        [0, -1, 0, 0, -1, 4, -1, 0, -1, 0],
        [0, 0, -1, 0, 0, -1, 4, 0, 0, 0],
        [0, 0, 0, 0, -1, 0, 0, 4, -1, -1],
        [0, 0, 0, 0, 0, -1, 0, -1, 4, 0],
        [0, 0, 0, 0, 0, 0, 0, -1, 0, 4]
    ])
    return A
    
    
def Assemble2dFEMmatrix(N):
    # Defining the matrices
    def T(M):
        return np.diag(-1 * np.ones(M - 1), -1) + np.diag(-1 * np.ones(M - 1), 1) + np.diag(4 * np.ones(M))

    def E(k, l):
        return np.eye(k, l)

    def O(D):
        return np.zeros((D, D))


    return A

