//**********************************************************************
// File Name: C5.1_31.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Compute the determinants of the first tel Hilbert matrices.
for i = 1:10
  d(i) = det(hilbertMatrix(i));
  end

// Here, we're interested in the pivots of hilbertMatrix(5);
[l,u] = lu(hilbertMatrix(5));
p = diag(u);

// Display results.
disp(d);
disp(p);





