//**********************************************************************
// File Name: Ci1.1-.sci
//**********************************************************************

exec('utils.sci',-1)

//**********************************************************************
// Mainline code.
//**********************************************************************ear
// Set display format to print to 3 decimal places. This is so FORTRAN!
format(10);

// Construct a 8 by 8 Hilbert matrix.
H = hilbertMatrix(8);

// Compute eigenvectors and eigenvalues.
[X,Lamda] = spec(H);

// Grab a of the diagonal entries of Lamda.
for k = 1:8
 lamdaKK(k) = Lamda(k,k);
end

// Force a column vectoe.
lamdaKK = lamdaKK(:);

// Compute minumim and maximum.
lamdaMin = min(lamdaKK);
lamdaMax = max(lamdaKK);

// Compute the norm of the inverse of H.
hInvNorm = norm(H \ 1);
