"""
MATH3036 CW1 DEtools module

@author: your name
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
    

    return b



def AssemblePGmatrixGenDE(N,c):
    """
    DESCRIPTION

    Parameters
    ----------
    N : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    
    Returns
    -------
    A: TYPE
       DESCRIPTION.

    """


    return A



def AssembleFEMmatrix(N,c):
    """
    DESCRIPTION

    Parameters
    ----------
    N : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    
    Returns
    -------
    A: TYPE
       DESCRIPTION.

    """
    

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

    
    
    return b






def AssembleMatrix_M(N):
    
    return M

    
def AssembleMatrix_K(N):
    
    return K


def AssembleVector_u0(N,u0_func):
        
    return u0_vec


def HeatEqFEM(tau,alpha,Ntime,N,u0_func):    
    u_array = np.zeros((Ntime+1, N-1))
    
    return u_array




def Assemble2dFEMvector(N):

    return b


def Assemble2dFEMmatrixN6():
        
    return A
    
    
def Assemble2dFEMmatrix(N):

    return A

