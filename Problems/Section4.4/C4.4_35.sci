//**********************************************************************
// File Name: C4.4_35.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct the initial matrix.
A = eye(4,4) - diag([1 1 1],-1);

// Perform the QR decomposition.
[Q,R] = qr(A,'e');

for k = 1:4
  // Compute normalization factors for the current column.
  n = abs(min(Q(:,k)));

  // Construct column of matrix with no normalization.
  O(:,k) = Q(:,k) / n;
end
