//**********************************************************************
// File Name: C3.4_41.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
p21 = [0 1 0; 1 0 0; 0 0 1];
p31 = [0 0 1; 0 1 0; 1 0 0];
p32 = [1 0 0; 0 0 1; 0 1 0];
p32_21 = p32 * p21;
p21_32 = p21 * p32;

// Fill in the identify elements plus the extras.
I = p21 + p31 + p32;

// Subtract out the extras.
I = I - p32_21 - p21_32;

// Take the transposes so that we can construct the coordinate vectors.
p21 = p21';
p31 = p31';
p32 = p32';
p32_21 = p32_21';
p21_32 = p21_32';

// Construct the coordinate row vectors.
r21 = p21(1:$)';
r31 = p31(1:$)';
r32 = p32(1:$)';
r32_21 = p32_21(1:$)';
r21_32 = p21_32(1:$)';

// Construct matrix of coordinate row vectors.
a = [r21; r31; r32; r32_21; r21_32];

// Compute the rank to determine the number of linearly independent rows.
r = rank(a);







