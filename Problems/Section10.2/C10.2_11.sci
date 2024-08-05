//**********************************************************************
// File Name: C10.2_11.sci
//**********************************************************************

//******************************************************************
// 
//  Name: FixedFixed
//
//  Purpose: The purpose of this function is to compute the
//  difference matrix that is used to generate the fixed-fixed
//  matrix, K = A'CA.
//
//  Calling Sequence: A = FixedFixed(n)
//
//  Inputs:
//
//    m - The number of columns in the output matrix.
//
//  Outputs:
//
//    A - the difference matrix.
//
//******************************************************************
function A = FixedFixed(n)

  // Construct identity matrix.
  A = eye(n,n);

  // Construct e -1 vector.
  o = ones(1,n-1) * -1;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct temporary matrix. The '-1' implies below the main
  //diagonal.
  D = diag(o,-1);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Add to our fixed-fixed difference matrixmatrix.
  A = A + D;
 
  // Append the final row.
  A(n+1,n) = -1;

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Crrate force vector.
f = ones(100,1) * 0.01;

// Construct fixed-fixed matrix.
a1 = FixedFixed(100);

// Construct fixed-free matrax.
a2 = a1;
a2($,:) = [];

// Create the K matrices.  C is an identity matrix so it is not used.
k1 = a1' * a1;
k2 = a2' * a2;

// Solve the two systems.
u1 = k1 \ f;
u2 = k2 \ f;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(211);
a = gca();
a.grid = [1 1];
title('Fixed-fixwd Solution');
plot(u1);

subplot(212);
a = gca();
a.grid = [1 1];
title('Fixed-free Solution');
plot(u2);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

