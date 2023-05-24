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
  // Extract indices of nonzero values in current column.
  idx = find(Q(:,k) <> 0);

  // Retrieve the nonzero values from in the current column.
  x = Q(idx,k);

  // Compute normalization factors for the current column.
  n = min(abs(x));

  // Construct column of matrix with no normalization.
  O(:,k) = Q(:,k) / n;
end
