//**********************************************************************
// File Name: C5.1_32.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set up the matrix size vector.
n = [50 100 200 400]';

// Compute the determinants.
for i = 1:4
  d(i) = det(rand(n(i),n(i)));
  d_n(i) = det(rand(n(i),n(i),'normal'));
end

// Display results.
disp([d d_n]);




