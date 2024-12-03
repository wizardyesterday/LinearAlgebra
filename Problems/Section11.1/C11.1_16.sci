//**********************************************************************
// File Name: C11.1_16.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create red-black permutation matrix.
P = permutationMatrix([1 3 5 7 9 2 4 6 8 10]);

// Create second difference matrix..
K = toeplitz([2 -1 zeros(1,8)]);

// Compute permuted red-black matrix.
A = P * K * P';

// Retrieve the upper right 5 by 5 block.
D = A(1:5,6:10);
