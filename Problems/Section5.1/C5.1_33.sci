//**********************************************************************
// File Name: C5.1_33.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************

// Compute the 5 by 5 matrix of 1's and -1's that has the maximum
// value of the determinant.
[h,d] = buildHadamardMatrix(5);

// Display results.
// Hadamard matrix.
disp(h);
// Determinant.
disp(d);

