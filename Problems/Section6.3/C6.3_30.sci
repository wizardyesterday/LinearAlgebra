//**********************************************************************
// File Name: C6.3_30.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set the step size.
a = %pi / 32;

// Construct the matrix.
A = [1-a^2  2*a; -2*a 1-a^2] / (1 + a^2);

// Set initial vector.
U = [1 0]';

for n = 1:33

  // Display current value.
  disp(U');

  // Update.
  U = A * U;

end

