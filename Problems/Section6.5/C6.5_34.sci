//**********************************************************************
// File Name: C6.5_34.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set parameter.
h = 1/6;

// Construct the second difference matrix.
S = [2 -1 0 0 0;
     -1 2 -1 0 0;
     0 -1 2 -1 0;
     0 0 -1 2 -1;
     0 0 0 -1 2];

// Contruct the sine and the eigenvalues of S*Q.
for k = 1:5
  // Populate the eigenvalue entry.
  lam(k) = 2 - (2 * cos(k * h * %pi));

  for i = 1:5
    // Populate the matrix entry.
    Q(i,k) = sin(i * k * h * %pi); 
  end
end

// Clean up the residual values.
Q = clean(Q);

// Compute the product.
A = clean(S * Q);

// Show that Q can be recovered by lamda.
q1 = A(:,1) / lam(1);
q2 = A(:,2) / lam(2);
q3 = A(:,3) / lam(3);
q4 = A(:,4) / lam(4);
q5 = A(:,5) / lam(5);

Qrecovered = [q1 q2 q3 q4 q5];






