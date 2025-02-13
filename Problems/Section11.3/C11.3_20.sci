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
//
//  Name: LanczosMethod
//
//  Purpose: The purpose of this function is to construct an n by n
//  orthogonal matrix using the Lanczos iteration method.
//
//
//  Calling Sequence: Q = LanczosMethod(A)
//
//  Inputs:
//
//    A - The input matrix.
//
//  Outputs:
//
//    Q - The constructed orthgonal matrix.
//
//**********************************************************************
function Q = LanczosMethod(A)

  // Create normalized vector.
  q1 = [1; -1; 0];
  q1 = norm(q1);

  // Set initial values.
  b = 1;
  r = q1;

  q2 = r / b;


  Q = [];

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create second difference matrix of order 3.
A = secondDifferenceMatrix(3);

// Create orthogonal matrix.
Q = LanczosMethod(A)



