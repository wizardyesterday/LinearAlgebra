//**********************************************************************
// File Name: C9.3-7.sci
//**********************************************************************

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct index vector.
t = 0:1;

// w for a 2-point Fourier transform..
w = exp(2*%pi*%i/2);

// Diagonal matrix for combining matrix.
D = diag(w.^t);

// 2-point Fourier matrix.
F2 = [1 1; 1 w];

// Consturt Combining matrix.
A = [eye(2,2) D; eye(2,2) -D];

// Construct block diagonal matrix.
B = [F2 zeros(2,2); zeros(2,2) F2];

// Construct permutation matrix.
P = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];

// Ffrst example matri.
c1 = [1 0 1 0]'

// second example matri.
c2 = [0 1 0 1]';

//The columns of this matrix contains the trajectory of c1.
S1(:,1) = c1;
S1(:,2) = P*c1;
S1(:,3) = B*P*c1;
S1(:,4) = A*B*P*c1;
S1 = clean(S1);

//The columns of this matrix contains the trajectory of c2.
S2(:,1) = c2;
S2(:,2) = P*c2;
S2(:,3) = B*P*c2;
S2(:,4) = A*B*P*c2;
S2 = clean(S2);








