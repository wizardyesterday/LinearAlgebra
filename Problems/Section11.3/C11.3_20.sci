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
  q1 = [1 -1 0]';
  q1 = q1 / norm(q1);

  // Construct initial orthogonal matrix. matix.
  Q = zeros(3,3);

  // Set initial column.
  Q(:,1) = q1;


endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create second difference matrix of order 3.
A = secondDifferenceMatrix(3);

// Create orthogonal matrix.
Q = LanczosMethod(A)



