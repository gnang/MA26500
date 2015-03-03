def gaussian_elimination(Cf, rs):
    """
    Outputs the row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: [RefA, c] = gaussian_elimination(Matrix(SR,2,2,[var('a'+str(i)+str(j)) for i in range(2) for j in range(2)] ), Matrix(SR,2,1,[var('b0'),var('b1')])
        sage: RefA
        [      1 a01/a00]
        [      0       1]
        sage: c
        [                                b00/a00]
        [(a10*b00/a00 - b10)/(a01*a10/a00 - a11)]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    A = copy(Cf); b = copy(rs)
    # Zero padding the matrix if necessary.
    if A.nrows() < A.ncols():
        # Temporary matrix for A
        Ta = Matrix(SR, zero_matrix(A.ncols(), A.ncols()))
        Ta[:A.nrows(), :A.ncols()] = copy(A)
        # Temporary matrix for b
        Tb=Matrix(SR, zero_matrix(A.ncols(), b.ncols())) 
        Tb[:b.nrows(), :b.ncols()] = copy(b)
        # replacing the matrix with the zero apdding.
        A=Ta; b=Tb
    # Initialization of the row and column index
    i=0; j=0;
    while i < A.nrows() and j < A.ncols():
        if (A[i:,j]).is_zero():
            # Incrementing the column index
            j=j+1
        if (A[i:,:].is_zero())==False:
            while (A[i,j]).is_zero(): 
                Ta=A[i:,:]
                Tb=b[i:,:]
                # Initializing the cyclic shift permutation matrix
                Id=identity_matrix(Ta.nrows())
                P=sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                Ta=P*Ta; Tb=P*Tb
                A[i:,:]=Ta
                b[i:,:]=Tb
            # Performing the row operations.
            b[i,:]=(1/A[i,j])*b[i,:]
            A[i,:]=(1/A[i,j])*A[i,:]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if A[r,:].is_zero():
                    r=r+1
                else:
                    b[r,:]=-A[r,j]*b[i,:]+b[r,:]
                    A[r,:]=-A[r,j]*A[i,:]+A[r,:]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A, b]

