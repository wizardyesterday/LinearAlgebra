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
//  Calling Sequence: [S,T] = GaussSeidelSplit(A)
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
//    T - The negative of the upper triangular portion of A.
//
//**********************************************************************
function [S,T] = GaussSeidelSplit(A)

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
//
//  Name: GaussSeidelIteration
//
//  Purpose: The purpose of this function is to perform the Gauss-Seidel
//  iteration to compute a solution to Ax = B.
//
//  Calling Sequence: [count,deltaXx,] =
//                        GaussSeidelIteration(S,T,b.tolerance)
//
//  Inputs:
//
//    S - The elements below and including the main diagonal of A.
//    elements.  This is the lower triangular portion of A, which was
//    previously split.
//
//    T - The negative of the upper triangular portion of A, which was
//    previously split.
//
//    b - The constand vector for Ax = b.
//
//  Outputs:
//
//    count - The number of iterations that occurred in order to reach a
//    solution.
//
//    deltaX - The norm of the difference of x between the final
//    iteration and the previous iteration.
//
//    deltaX - The norm of the difference of x between the final
//    iteration and the previous iteration.
//
//    x - The estimated solution.
//
//**********************************************************************
function [count,x,deltaX] = GaussSeidelIteration(S,T,b,tolerance)

  // Ensure we have a column vector.
  b = b(:);

  // Save the length of b for creation of x.
  len = length(b);

  // Compute iterator mastrix and inverse matrix.
  B = S \ T;
  Sinverse = S \ 1;

  // Create initial guess for solution.
  x = zeros(len,1);

  // Set up for loop wntry.
  done = 0;
  count = 0;

  while done == 0 then

    // One more iteration has occurred.
    count =  count + 1;

    // Compute estimate.
    xNew = (B * x) + (Sinverse * b); 

    //Compute norm of error.
    deltaX = xNew - x;
    deltaX = norm(deltaX);

    // Update estimate.
    x = xNew;

    if deltaX <= tolerance then
      // Bail out of loop, we're done.
      done = 1;
    else
      if count > 1000
        // Bail out of loop, we failed to converge.
        done = 1;
      end
    end

  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Compute right-side vectors.
b10 = zeros(10,1);
b10(1) = 1;
b20 = zeros(20,1);
b10(1) = 1;
b10 = zeros(50,1);
b50(1) = 1;

// Create second difference matrices of order 10, 20, 50..
A10 = toeplitz([2 -1 zeros(1,8)]);
A20 = toeplitz([2 -1 zeros(1,18)]);
A50 = toeplitz([2 -1 zeros(1,48)]);


// Create the Splits.
[S10,T10] = GaussSeidelSplit(A10);
[S20,T20] = GaussSeidelSplit(A20);
[S50,T50] = GaussSeidelSplit(A50);

b = zeros(2,1);
b(1) = 0;
a = [1 2; 6 7];
[s,t] = GaussSeidelSplit(a);
//function [count,x,deltaX] = GaussSeidelIteration(S,T,b,tolerance)


