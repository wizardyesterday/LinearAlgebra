//**********************************************************************
// File Name: C2.2_34.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Generate random vector of length 3 with unit variance.
v = generateGaussianProcess(1,3,1);

// Force column vector.
v = v(:);

// Construct a unit vector.
u = v / norm(v);

// Generate 30 realizations of random vectors with unit variance.
V = generateGaussianProcess(30,3,1);

// Construct unit column vectors.
for j = 1:30
  U(:,j) = V(:,j) / norm(V(:,j));
end

// Compute vector of dot products.
uDotU = u' * U;

// Compute the average of |uU|_j, j = 1, ... ,30.
uDotU_avg = mean(abs(uDotU));

// Output result.
printf("Mean value of |u dot U|: %f\n",uDotU_avg);







