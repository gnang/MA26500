# Implementation of the Gaussian Elimination
# The function expect as input a symbolic matrix A
def Elimination(A):
    # Matrix size parameters.
    m=A.nrows(); n=A.ncols()
    # Carrying the first part of the Gaussian elimination steps
    for i in range(1,m):
        A[m-i,:] = sum([exp(I*2*pi*k/(m-i+1))*prod([A[j,0] for j in range(m-i+1) if j!=k])*A[k,:] for k in range(m-i+1)])
    # For displying purposes 
    return A

def GaussianElimination(A):
    # Matrix size parameters.
    m=A.nrows(); n=A.ncols()
    for i in range(A.nrows()-1):
        A[i:,i:]=Elimination(A[i:,i:])
    return A
    
