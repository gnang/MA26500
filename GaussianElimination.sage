# Short implementation of Gaussian Elimination
# perform row operation on the matrix to create
# zeros bellow the pivot locations
def Elimination(M):
    A = copy(M)
    # Matrix size parameters.
    m=A.nrows(); n=A.ncols()
    # Performing row operations
    # to create zeros bello the
    # position.
    for i in range(1,m):
        A[m-i,:] = sum([exp(I*2*pi*k/(m-i+1))*prod([A[j,0] for j in range(m-i+1) if j!=k])*A[k,:] for k in range(m-i+1)])
    return Matrix(SR, m, n, [(A[i,j]).simplify_full() for i in range(m) for j in range(n)])

# Computes an almost row echelon form by avoinding divisions
def GaussianElimination(A):
    # Matrix size parameters.
    m=A.nrows(); n=A.ncols()
    # Performs Gaussian elimination by recursively calling
    # the elimination function implemented above.
    for i in range(A.nrows()-1):
        A[i:,i:]=Elimination(A[i:,i:])
    return A

def SymbolicLinearSolver(A,b):
    # Initializing the size parameter
    sz = b.nrows()
    # Initializing the variables
    x=Matrix(SR,sz,1,[var('x'+str(i)) for i in range(sz)])
    # Initializing the matrix on which we will permform gaussian elimination.
    M=Matrix(SR,zero_matrix(sz,sz+1)); M[0:sz,0:sz]=A[0:sz,0:sz]; M[:,sz]=b[:,0]
    # Computing the almost row echelon form.
    Rslt=GaussianElimination(M); v=Rslt[0:sz,0:sz]*x+Rslt[:,sz]
    # Deriving the solution
    Solutions = [[] for i in range(sz)]
    Solutions[sz-1].append(var('x'+str(sz-1)))
    Solutions[sz-1].append((((v[sz-1,0])/(v[sz-1,0]).coeff(var('x'+str(sz-1)))).simplify_full()).subs( dict([(var('x'+str(sz-1)),0)])  ))
    for i in range(1,sz):
        Solutions[sz-i-1].append(var('x'+str(sz-i-1)))
        f = v[sz-i-1,0].subs(dict([(var('x'+str(j)), -Solutions[j][1]) for j in range(sz-i,sz)]))
        Solutions[sz-i-1].append((f/f.coeff(var('x'+str(sz-i-1)))).subs( dict([(var('x'+str(sz-i-1)),0)]) ))
    return Solutions
 
