//**********************************************************************
// File Name: C7.1_9.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct uniformly distributed random matrix.
A = rand(20,20,'uniform');

// Construct normally distributed random matrix.
B = rand(20,20,'normal');

subplot(211);
title('Singular Values for Uniform Distribution');
plot(svd(A));

subplot(212);
title('Singular Values for Normal Distribution');
plot(svd(B));

