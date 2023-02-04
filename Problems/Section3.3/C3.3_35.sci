//**********************************************************************
// File Name: C3.3_35.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: secondDifferenceMatrix
//
//  Purpose: The purpose of this function is to construct a second
//  difference matrix.  A second difference matrix is an n by n (square)
//  matrix with 2's along the main diagonal, and -1's above and below the
//  main diagonal.
//
//  Calling Sequence: A = secondDifferenceMatrix(n)
//
//  Inputs:
//
//    n - The number of rows and columns of the second difference
//    matrix.
//
//  Outputs:
//
//    A - An n by n second difference matrix.
//
//**********************************************************************
function A = secondDifferenceMatrix(n)

  // Construct identify entries of the matrix.
  A = 2 * eye(n,n);

  for i = 2:n
    // Fill in entries above the main diagonal.
    A(i-1,i) = -1;

    // Fill in entries below the main diagonal.
    A(i,i-1) = -1;
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct a 9 by 9 second difference matrix.
K = secondDifferenceMatrix(9);

// Construct a 9 by 1 matrix of 10's.
b = 10 * ones(9,1);

// Solve Ax = b, for x.
x = K \ b;



