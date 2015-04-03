def Gram_Schmidt(M):
    # Initialization of the symbolic matrix to be used as output
    Rs = Matrix(SR, zero_matrix(M.nrows(), M.ncols()))

    # This implementation assumes that the first column is a non zero column.
    # Initializing the first vector to a unit vector.
    Rs[:,0] = (1/sqrt(sum([M[j,0]^2 for j in range(M.nrows())]))) * M[:,0]

    # Loop removing the bad components for all the remaining vectors
    for i in range(1, M.ncols()):
        v = M[:,i]
        # Removing the bad components for the column M[:,i] 
        v = v - sum([(v.transpose()*Rs[:,j])[0,0]*Rs[:,j] for j in range(i)])

        # Checking that the resulting column is not the zero column
        if v != Matrix(SR,zero_matrix(1,M.nrows())):
            Rs[:,i] = (1/sqrt(sum([v[j,0]^2 for j in range(M.nrows())]))) * v
        else:
            Rs[:,i] = v
    return Rs
