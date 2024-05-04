//**********************************************************************
// File Name: C9.3-8.sci
//**********************************************************************

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct index vector.
t = 0:3;

// w for a 4-point Fourier transform..
w = exp(2*%pi*%i/4);

// Diagonal matrix for combining matrix.
D = diag(w.^t);

// 4-point Fourier matrix.
F4 = [1 1 1 1; 1 w w^2 w^3; 1 w^2 w^4 w^6; 1 w^3 w^6 w^9];

// Consturt Combining matrix.
A = [eye(4,4) D; eye(4,4) -D];

// Construct block diagonal matrix.
B = [F4 zeros(4,4); zeros(4,4) F4];

// Construct permutation matrix.
P = [1 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 1 0;
     0 1 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 1];

// Ffrst example matri.
c1 = [1 0 1 0 1 0 1 0]'

// second example matri.
c2 = [0 1 0 1 0 1 0 1]';

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








