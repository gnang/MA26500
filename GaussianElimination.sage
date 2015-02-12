def gaussian_elimination(Cf, rs):
    A=copy(Cf); b=copy(rs)
    # Zero padding the matrix if necessary.
    if A.nrows() < A.ncols():
        # Temporary matrix for A
        Ta=Matrix(SR, zero_matrix(A.ncols(), A.ncols()))
        Ta[:A.nrows(), :A.ncols()] = copy(A)
        # Temporary matrix for b
        Tb=Matrix(SR, zero_matrix(A.ncols(), b.ncols())) 
        Tb[:b.nrows(), :b.ncols()] = copy(b)
        # replacing the matrix with the zero padding.
        A = Ta; b = Tb
    # Initializing the cyclic shift permutation matrix
    Id = identity_matrix(A.nrows())
    P  = sum([Id[:,k]*Id[mod(k+1,A.nrows()),:] for k in range(A.nrows())])
    # Initialization of the row and column index
    i = 0; j = 0;
    while i < A.nrows() and j < A.ncols():
        if (A[:,j]).is_zero():
            # Incrementing the column index
            j = j + 1
        else:
            while (A[i,j]).is_zero():
                A = P*A; b = P*b
            # Performing the row operations.
            b[i,:] = (1/A[i,j])*b[i,:]
            A[i,:] = (1/A[i,j])*A[i,:]
            for r in range(i+1,A.nrows()):
                b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                A[r,:] = -A[r,j]*A[i,:]+A[r,:]
        # Incrementing the row and column index.
        i = i + 1; j = j + 1
    return [A, b]

