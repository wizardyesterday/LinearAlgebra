//**********************************************************************
// File Name: Ci1.1-.sci
//**********************************************************************

exec('utils.sci',-1)

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct b with ||b+ = 1.
b = ones(9,1) / 3;

// Construct a 9 by 9 Hilbert matrix.
H = hilbertMatrix(9);

// Compute eigenvectors and eigenvalues.
[X,Lamda] = spec(H);

// Grab a of the diagonal entries of Lamda.
for k = 1:9
 lamdaKK(k) = Lamda(k,k);
end

// Force a column vectoe.
lamdaKK = lamdaKK(:);

// Compute minumim and maximum.
lamdaMin = min(lamdaKK);
lamdaMax = max(lamdaKK);

//Compute the norms of eaxh column of X.
for k = 1:9
  HNorm(:,k) = norm(H(:,k));
end

// Force a column vectoe.
hNorm = HNorm(:);

// Compute the solution to Hx = b.
x = H \ b;

// Compute ||x||.
xNorm = norm(x);




