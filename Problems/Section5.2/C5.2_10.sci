//**********************************************************************
// File Name: C5.2_10.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Compute the permutations of (1,2,3,4).
p = perms(1:4);

// Construct the 4 by 4 identify matrix.
I4 = eye(4,4);

// Iniitalize to the first element of the output vector.
evenResultIndex = 1;

// Start with a clean slate.
result = 0;

for i = 1:24
  // Compute a candidate permutation matrix.
  M = permutationMatrix(p(i,:));

  // Compute the determinant of the permutation matrix.
  d = det(M);

  if d == 1
    // We have an even permutation.
    result(evenResultIndex) = det(I4 + M);

    // Reference the next element to be stored.
    evenResultIndex = evenResultIndex + 1;
  end
end

// Display results.
disp(result);
