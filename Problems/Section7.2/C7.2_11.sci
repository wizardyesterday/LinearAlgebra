//**********************************************************************
// File Name: C7.2_11.sci
//**********************************************************************

//**********************************************************************
// Mainline code.
//**********************************************************************
// Conwtruct origin vector.
o = [0 0];

// Conwtruct some arbitrary unit vector.
u = [1 1] / sqrt(2);

// Construct problem matrix.
A = [1 1;1 0];

// Retrieve the eigenvalues and eigenvectors of A.
[X,D] = spec(A);

// We're only interested in the first eigenvector.
x1 = X(:,1);

// Construct the data for the first plots.
z1 = [o' x1];
z2 = [o' A*x1];

// Construct the data for the second plots.
z3 = [o' u'];
z4 = [o' A*u'];

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Perform eigenvector plot.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(211);
title('Plot of Eigenvector x1 and Ax1');
a = gca();
a.grid = [1 1];

plot(z1(1,:), z1(2,:),"r");
plot(z2(1,:), z2(2,:),"b");

legend("x1","Ax1");
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Perform arbitrary unit vector plot.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(212);
title('Plot of Arbitrary Unit Vector u and Au');

a = gca();
a.grid = [1 1];

plot(z3(1,:), z3(2,:),"r");
plot(z4(1,:), z4(2,:),"b");

legend("u","Au");
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

