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

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct initial orthogonal matrix. with initial column vector.
  // Q contains column vectors: [q1 q2 ... qn'.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Preallocate the matrix.
  Q = zeros(n,n);

  // q1 = [1 0 ... 0]'.
  Q(1) = 1;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Set initial column, r0 equal to q1.
  r = Q(:,1);

  // Chose b0 equal to 1.
  b(1) = 1;

   // set initial conditions for loop.
   done = 0;
   j = 1;

  while done == 0

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // A single Lanczos iteration.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // q_j+1 = r_j / b_j.
    Q(:,j+1) = r / b(j);

    // Increment.
    j = j + 1;

    // a_j = q_j' * A * q_j.
    a = Q(:,j)' * A * Q(:,j);

    // r_j = A * q_j - b_j-1 * q_j-1 - a_j * q_j.
    r = A * Q(:,j) - b(j-1) * Q(:,j-1) - a * Q(:,j);

    // b_j = ||r_j||.
    b(j) = norm(r);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

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



