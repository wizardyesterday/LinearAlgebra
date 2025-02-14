//**********************************************************************
// File Name: C11.3_20.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Functions.
//**********************************************************************
//
//**********************************************************************
// Mainline code.
//**********************************************************************
// Create second difference matrix of order 8.
A = secondDifferenceMatrix(8);

// Create orthogonal matrix.
Q = LanczosMethod(A);

// Display result.  It should be an identity matrix.
disp(clean(Q'*Q));



