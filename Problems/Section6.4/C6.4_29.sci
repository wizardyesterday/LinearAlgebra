//**********************************************************************
// File Name: C6.4_29.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set the step size.
t = -08:0.01:8;

// Set up matrices.
A = [1 0; 0 2];
B = [8 1; 1 0];

for i = 1:length(t)
  // Compute the function of t.
  C = A + (t(i) * B);

  // Compute eigenvectors.
  L = spec(C);

  // Separate the eigenvalue.
  L1(i) = L(1);
  L2(i) = L(2);
end

// Display results.
disp([t' L1 L2 (L1 - L2)]);

// Plot results.
plot(t,(L1 - L2));
title('Eigenvalue Difference (Lamda1 - Lamda2)');

