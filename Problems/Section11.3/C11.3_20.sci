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
//  Calling Sequence: Q = LanczosMethod(n)
//
//  Inputs:
//
//    n - The order of the orthogonal matrix.
//
//  Outputs:
//
//    Q - The constructed orthgonal matrix.
//
//**********************************************************************
function Q = LanczosMethod(n)

  // Create normalized vector.
  q1 = [1; -1; 0];
    q1` = norm(q1);

  // Set initial values.
  b = 1;
  r = q1;

  Q = [];

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create second difference matrix of order 3.
A = secondDifferenceMatrix(3);

// Create orthogonal matrix.
Q = LanczosMethod(3)



