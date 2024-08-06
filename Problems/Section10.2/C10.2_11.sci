//**********************************************************************
// File Name: C10.2_11.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Crrate force vector.
f = ones(100,1) * 0.01;

// Construct fixed-fixed matrix.
a1 = FixedFixed(100);

// Construct fixed-free matrax.
a2 = a1;
a2($,:) = [];

// Create the K matrices.  C is an identity matrix so it is not used.
k1 = a1' * a1;
k2 = a2' * a2;

// Solve the two systems.
u1 = k1 \ f;
u2 = k2 \ f;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(121);
a = gca();
a.grid = [1 1];
title('Fixed-fixwd Solution');
xlabel('mass number');
ylabel('movement');
plot(u1);

subplot(122);
a = gca();
a.grid = [1 1];
title('Fixed-free Solution');
xlabel('mass number');
ylabel('movement');
plot(u2);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

