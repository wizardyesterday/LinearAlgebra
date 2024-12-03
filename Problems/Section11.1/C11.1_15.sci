//**********************************************************************
// File Name: C11.1_15.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create Toeplitz matrix.
K = toeplitz([2 -1 zeros(1,8)]);

// Permute eows randomly.
KK = K(randperm(10),randperm(10));

// Perform LU factorizations
[L,U] = lu(K);
[LL,UU] = lu(KK);

// Count zeros.
zL = nnz(L);
zLL = nnz(LL);
