//**********************************************************************
// File Name: C2.2_29.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Initialize accumulators.
p = zeros(1,3);

for j = 1:100
  // Compute LU decomposition
  [L,U] = lu(rand(eye(3,3)));

  for k = 1:3
    // Update accumulators.
    p(k) = U(k,k);
  end
end

// Compute mean values of the absolute values of the pivots.
p = p / 100;

// Output results.
for k = 1:3
  printf("|p%d%d|: %f\n",k,k,abs(p(k)));
end



