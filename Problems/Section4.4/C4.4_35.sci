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

// Compute normalization factors for each column.
n1 = abs(min(Q(:,1)));
n2 = abs(min(Q(:,2)));
n3 = abs(min(Q(:,3)));
n4 = abs(min(Q(:,4)));

// Compute columns with no normalization.
A1 = Q(:,1) / n1;
A2 = Q(:,2) / n2;
A3 = Q(:,3) / n3;
A4 = Q(:,4) / n4;

// Construct the final matrix with orthogonal columns (not orthonormal).
O(:,1) = Q(:,1) / n1;
O(:,2) = Q(:,2) / n2;
O(:,3) = Q(:,3) / n3;
O(:,4) = Q(:,4) / n4;
