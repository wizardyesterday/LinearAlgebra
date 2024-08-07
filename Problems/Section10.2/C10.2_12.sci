//**********************************************************************
// File Name: C10.2_12.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Crrate force vector.
f = ones(7,1);

// Construct matrix with ones above the main diagonal.
E = diag(ones(6,1),1);

// Construce our difference matrix.
K = 2 * eye(7,7);
K = K - E - E';
K = K * 64;

// Construct another difference matrix.
D = E - eye(7,7);
D = D * 80;
// Compute forward diference solution.
uForward = (K + D) \ f;

// Compute backward difference solution.
uBackward = (K - D') \ f;

// Compute centered difference solution.
uCentered = (K + D/2 - D'/2) \ f;


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(131);
a = gca();
a.grid = [1 1];
title('Forward Difference Solution.');
plot(uForward);

subplot(132);
a = gca();
a.grid = [1 1];
title('Backward Difference Solution.');
plot(uBackward);

subplot(133);
a = gca();
a.grid = [1 1];
title('Centered Difference Solution.');
plot(uCentered);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

