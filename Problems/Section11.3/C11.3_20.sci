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

  // Compute order of A.
  n = size(A);
  n = n(1);

  // Construct initial orthogonal matrix. with initial column vector.
  Q = zeros(n,n);
  Q(1) = 1;

  // Set initial column, r = q1.
  r = Q(:,1);
  b(1) = 1;

   // set initial conditions for loop.
   done = 0;
   j = 1;

   while done == 0

    Q(:,j+1) = r / b(j);
    j = j + 1;
    a = Q(:,j)' * A * Q(:,j);
    r = A * Q(:,j) - b(j-1) * Q(:,j-1) - a * Q(:,j);
    b(j) = norm(r);

    if j == n then
      // bail out.
      done = 1;
    end

  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create second difference matrix of order 3.
A = secondDifferenceMatrix(8);

// Create orthogonal matrix.
Q = LanczosMethod(A)



