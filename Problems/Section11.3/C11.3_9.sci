//**********************************************************************
// File Name: C11.3_9.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Functions.
//**********************************************************************
//
//**********************************************************************
//
//  Name: GaussSeidelSplit
//
//  Purpose: The purpose of this function is to split a square matrix
//  using the Gauss-Seidel factorization.
//
//
//  Calling Sequence: [ST] = GaussSeidelSplit(A)
//
//  Inputs:
//
//    A - The matrix to be split.
//
//  Outputs:
//
//    S - The elements below and including the main diagonal of A.
//    elements.  This is the lower triangular portion of A.
//
//  T - The negative of the upper triangular portion of A.
//
//**********************************************************************
function [S, T] = GaussSeidelSplit(A)

  // Compute the order of A.
  N = size(A);
  N = N(1);

  // Preallocate result matrices.
  S = zeros(N,N);
  T = zeros(N,N);

  // Retrieve lower triangular portion of A.
  for i = 1:N
    for j = i:N
      S(j,i) = A(j,i);
    end
  end

  // Compute -(upper triangular portion of A).
  T = S - A;

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create 10 by 10 second difference matrix.
A10 = toeplitz([2 -1 zeros(1,8)]);

// Create the Splits.
[S10,T10] = GaussSeidelSplit(A10);

// Create 20 by 20 second difference matrix.
A20 = toeplitz([2 -1 zeros(1,18)]);

// Create the Splits.
[S20,T20] = GaussSeidelSplit(A20);

// Create 50 by 50 second difference matrix.
A50 = toeplitz([2 -1 zeros(1,48)]);

// Create the Splits.
[S50,T50] = GaussSeidelSplit(A50);

a = [1 2 3; 4 5 6; 7 8 9]
[S,T] = GaussSeidelSplit(a);

