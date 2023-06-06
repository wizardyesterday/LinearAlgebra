//**********************************************************************
// File Name: C5.5_31.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct the first ten Hilbert matrices.
h1 = hilbertMatrix(1);
h2 = hilbertMatrix(2);
h3 = hilbertMatrix(3);
h4 = hilbertMatrix(4);
h5 = hilbertMatrix(5);
h6 = hilbertMatrix(6);
h7 = hilbertMatrix(7);
h8 = hilbertMatrix(9);
h9 = hilbertMatrix(9);
h10 = hilbertMatrix(10);

// Compute the determinants of the Hilbert matrices.
d(1) = det(h1);
d(2) = det(h2);
d(3) = det(h3);
d(4) = det(h4);
d(5) = det(h5);
d(6) = det(h6);
d(7) = det(h7);
d(8) = det(h9);
d(9) = det(h9);
d(10) = det(h10);

// Here, we're interested in the pivots of hilbertMatrix(5);
[l,u] = lu(h5);
p = diag(u);
oneOverP = 1 ./ p;

// Display results.
disp(clean(d));
disp(p);
disp(oneOverP);





