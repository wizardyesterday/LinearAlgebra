//******************************************************************
// cab  A = c a b echelon factorization.
//
// [c, a, b] = cab(A) gives echelon bases for the column space in c
// and the row space in b
// b contains the nonzero rows of the echelon form rref(A)
// c contains the nonzero columns of the echelon form rref(A')'
// All extra nonzeros are below I in c and to the right of I in b
// a is the nonsingular submatrix formed by the pivot columns and 
// pivot rows of A.  Those columns of b and rows of c contain I.
//
// See also elim, rref.
//******************************************************************
function [c, a, b] = cab(A)

  [R, pivcol] = rref(A);
  [S, pivrow] = rref(A');
  b = R(1:rank(A), : );
  c = S(1:rank(A), : )';
  a = A(pivrow, pivcol);

endfunction

//******************************************************************
// cofactor  Matrix of cofactors.
//
// C = cofactor(A) returns the matrix of cofactors of A.
// If A is invertible, then inv(A) = C' / det(A).
//
// C = cofactor(A, i, j) returns the cofactor of 
// row i, column j of A.
//******************************************************************
function C = cofactor(A, i, j)

  if argn(2) == 3
  // Remove row i and column j to produce the minor.
    M = A;
    M(i,:) = [];
    M(:,j) = [];
    C = (-1)^(i+j) * det(M);
  else
    [n,n] = size(A);
    for i = 1:n
      for j = 1:n
     C(i,j) = cofactor(A, i, j);
      end
   end
  end

endfunction

//******************************************************************
// colbasis  Basis for the column space. 
//
// C = colbasis(A) returns the r pivot columns of A
// as a basis for the column space of A.
//
// See also fourbase.
//******************************************************************
function C = colbasis(A)

  [R, pivcol] = rref(A);
  C = A(:, pivcol);

endfunction

//******************************************************************
// cramer  Solve the system Ax=b.
// The matrix A is square and invertible.
//
// x = cramer(A, b) solves the square system Ax = b.
//******************************************************************
function x = cramer(A, b)

  [m, n] = size(A);
  if m ~= n
    error('Matrix is not square.') 
  end
  if det(A) == 0
    error('Matrix is singular.')
  end
  for j = 1:n
    B = A;
    B(:, j) = b;
    x(j) = det(B) / det(A);
  end
  x = x';

endfunction

//******************************************************************
// determ  Matrix determinant from plu.
//
// det = determ(A) computes the determinant of the square matrix A
// as the sign of the permutation times the product of pivots.
//******************************************************************
function det = determ(A)

  [P, L, U, sign] = splu(A);
  det = sign * prod(diag(U));

endfunction

//******************************************************************
// eigen2  Characteristic polynomial, eigenvalues, eigenvectors 
// of a 2 by 2 matrix.
//
// eigen2(A) prints the characteristic polynomial det(A-e*I),
// eigenvalues, and eigenvectors of A. 
//
// If A is not diagonalizable, its single eigenvector is 
// printed twice.
//******************************************************************
function eigen2(A)

  d = A(1,1)*A(2,2) - A(1,2)*A(2,1);
  t = A(1,1) + A(2,2);
  e1 = (t + sqrt(t^2 - 4*d))/2;
  e2 = (t - sqrt(t^2 - 4*d))/2;
  if A(1,2) ~= 0
     x1 = [A(1,2); e1-A(1,1)];
     x2 = [A(1,2); e2-A(1,1)];
  elseif A(2,1) ~= 0
     x1 = [e1-A(2,2); A(2,1)];
     x2 = [e2-A(2,2); A(2,1)];
  else
   x1 = [1; 0];
   x2 = [0; 1];
  end

  disp(' ')
  disp('For this matrix, the polynomial whose roots are the eigenvalues is:')
  disp(['   e^2 - ' num2str(t) '*e + ' num2str(d) ' = 0'])

  disp(' ')
  disp('The first eigenvalue and eigenvector are:')
  e1
  x1

  disp(' ')
  disp('The second eigenvalue and eigenvector are:')
  e2
  x2

endfunction

//******************************************************************
// eigval  Eigenvalues and their algebraic multiplicity.
//
// evalues = eigval(A) returns the distinct eigenvalues of A,
// with duplicates removed and sorted in decreasing order.
//
// [evalues, repeats] = eigval(A) also returns the row vector
// repeats that gives the multiplicity of each eigenvalue.
// The sum of the multiplicities is n.
//
// Examples: Let A = eye(n) and B = diag([3 4]).
// For A, evalues is 1 and repeats is n.
// For B, evalues is [4; 3]  and repeats is [1 1].
//******************************************************************
function [evalues, repeats] = eigval(A)

  tol = sqrt(%eps);
  lambda = sort(spec(A));
  lambda = round(lambda/tol) * tol;
  //
  // lambda gives all n eigenvalues (repetitions included).
  //
  evalues = unique(lambda);
  evalues = flipud(evalues);
  n = length(lambda);
  d = length(evalues);
  A = ones(n, 1) * evalues';
  B = lambda * ones(1, d);
  MATCH = abs(A-B) <= tol;
  //
  // MATCH is an n by d zero matrix except
  // MATCH(i,j) = 1 when lambda(i) = evalues(j).
  // Summing the columns gives the row vector repeats.
  //
  repeats = sum(MATCH);

endfunction

//******************************************************************
// eigvec  Eigenvectors and their geometric multiplicity.
//
// S = eigvec(A) returns the largest possible set of linearly
// independent eigenvectors of A. 
//
// [S, D] = eigvec(A) also returns the corresponding eigenvalues
// in the diagonal matrix D.
// Each eigenvalue in D is repeated according to the number of its
// linearly independent eigenvectors. This is its geometric multiplicity.
//
// Always A*S = S*D. If S is square then A is diagonalizable and
// inv(S)*A*S = D = LAMBDA.
//******************************************************************
function [S, D] = eigvec(A)

  [m, n] = size(A);
  I = eye(n,n);
  [evalues, repeats] = eigval(A);
  S = []; d = [];
  for k = 1 : length(evalues);
    s = nulbasis(A - evalues(k)*I);
    [ms, ns] = size(s);
    S = [S s];
    temp = ones(ns, 1) * evalues(k);
    d = [d; temp];
  end
  D = diag(d);

endfunction

//******************************************************************
// elim  E*A = R factorization.
//
// E = elim(A) returns the elimination matrix E
// that gives the reduced row echelon form E*A = R.
// If A is square and invertible, then E = inv(A).
//
// [E, R] = elim(A) returns the elimination matrix E 
// and the reduced row echelon form R.
//
// See also lu, slu, splu, plu.
//******************************************************************
function [E, R] = elim(A)

  [m, n] = size(A);
  I = eye(m,m);
  //
  // Elimination on the augmented matrix [A I] yields [R E].
  //
  RE = rref([A I]);
  R = RE(:, 1:n);
  E = RE(:, (n+1):(m+n));

endfunction

//******************************************************************
// findpiv  Used by plu to find a pivot for Gaussian elimination.
//
// [r, p] = findpiv(A(k:m, p:n), k, p, tol) finds the first element in
// the specified submatrix which is larger than tol in absolute value.
// It returns indices r and p so that A(r, p) is the pivot.
//******************************************************************
function [k, p] = findpiv(A, k, p, tol)

  [m,  n] = size(A);
  r = find(abs(A(:)) > tol);
  if isempty(r), return, end
  //
  r = r(1);
  j = fix((r-1)/m)+1;
  p = p+j-1;
  k = k+r-(j-1)*m-1;

endfunction

//******************************************************************
// grams  Gram-Schmidt orthogonalization of the columns of A.
// The columns of A are assumed to be linearly independent.
//
// Q = grams(A) returns an m by n matrix Q whose columns are 
// an orthonormal basis for the column space of A.
//
// [Q, R] = grams(A) returns a matrix Q with orthonormal columns
// and an invertible upper triangular matrix R so that A = Q*R.
//
// Warning: For a more stable algorithm, use [Q, R] = qr(A, 0) .
//******************************************************************
function [Q, R] = grams(A)

  [m, n] = size(A);
  Asave = A;
  for j = 1:n
    for k = 1:j-1
      mult = (A(:, j)'*A(:, k)) / (A(:, k)'*A(:, k));
      A(:, j) = A(:, j) - mult*A(:, k);
    end
  end
  for j = 1:n
    if norm(A(:, j)) < sqrt(eps)
      error('Columns of A are linearly dependent.')
    end
    Q(:, j) = A(:, j) / norm(A(:, j));
  end
  R = Q'*Asave;
endfunction

//******************************************************************
// X = house stores the "house" data set in X.
//
// The "house" data set is for use with plot2d.
// Try plot2d(A*X) for various 2 by 2 matrices A.
//
// See also hand.
//******************************************************************
function X = house

  X = [ -6  -6  -7   0   7   6   6  -3  -3   0   0  -6
        -7   2   1   8   1   2  -7  -7  -2  -2  -7  -7 ];

endfunction

//******************************************************************
// inverse  Matrix inverse by Gauss-Jordan elimination.
//
// Ainv = inverse(A) computes the inverse of A, if it exists.
//
// Row reduction applied to [A I] using elim produces [I Ainv].
//
// See also inv, elim. 
//******************************************************************
function Ainv = inverse(A)

  [m, n] = size(A);
  r = rank(A);
  if (r == m) & (r == n) 
    [Ainv, R] = elim(A);
  else
    Ainv = [];
    disp('Warning: A is not a square, invertible matrix.');
  end

endfunction

//******************************************************************
// leftnull  Basis for the left nullspace.
//
// LN = leftnull(A) returns a basis for the 
// left nullspace in the *columns* of LN.
//
// The left nullspace of A is the nullspace of A'.
// The command fourbase(A) finds a different basis
// for the left nullspace of A. 
//
// See also fourbase.
//******************************************************************
function LN = leftnull(A)

  LN = nulbasis(A');

endfunction

//******************************************************************
// linefit  Plot the least squares fit by a line.
//
// linefit(t, b), where t and b are vectors of the same length,
// displays the best line fit to the data points (t(i), b(i)).
//******************************************************************
function linefit(t, b)

  // We now insure that t and b are column vectors.
  t = t(:); b = b(:);

  // Form the matrix whose first column is all ones and
  // whose second column is the vector t.
  n = length(t);
  e = ones(n, 1);
  A = [e t];

  // Solve the least squares problem, A*x ~= b.
  xhat = lsq(A, b);
  c = xhat(1);
  d = xhat(2);

  // Plot the results.
  tline = [1.1*min(t)-0.1*max(t), 1.1*max(t)-0.1*min(t)];
  yline = c + d*tline;
  plot(t,b,'ro',t,c+d*t,'k*',tline,yline,'k-')
  if d >= 0 , sign = ' + '; else, sign = ' - '; end
  title(['Best line is ' num2str(c) sign num2str(abs(d)) '*t.'])
  xlabel('t')

endfunction

//******************************************************************
// lsq  Least squares solution to Ax=b.
// The columns of A are linearly independent.
//
// [xhat, p, e] = lsq(A, b) finds a least squares
// solution to the overdetermined system Ax ~= b.
// xhat solves the normal equations A'*A*xhat = A'*b.
// p is the projection of b onto the column space.
// e = b - p.
//******************************************************************
function [xhat, p, e] = lsq(A, b)

  if det(A'*A) == 0
    error('Columns of A are linearly dependent.')
  end
  xhat = partic(A'*A, A'*b);
  p = A * xhat;
  e = b - p;

endfunction

//******************************************************************
// normal  Eigenvalues and eigenvectors of a normal matrix A.
//
// U = normal(A) returns a set of orthonormal eigenvectors for A.
//
// [U, LAMBDA] = normal(A) also returns the corresponding eigenvalues 
// on the diagonal of LAMBDA. The eigenvalues on the diagonal of 
// LAMBDA are sorted by magnitude.
//
// Normal matrices (A'*A = A*A') may have complex eigenvalues and 
// eigenvectors. If A itself is complex, A' is its conjugate transpose.
//
// See also eig, eigval, eigvec, symmeig.
//******************************************************************
function [U, LAMBDA] = normal(A)

  [m, n] = size(A);
  E = A'*A - A*A';
  if norm(E) <= sqrt(%eps)
  //
  // The eigenvectors in S are linearly independent but not orthogonal.
  // Eigenvectors for different eigenvalues *are* orthogonal.
  // Gram-Schmidt (qr(S)) gives orthonormal eigenvectors in U.
  //
    [S, LAMBDA] = eigvec(A);
    [U, R] = qr(S);
  else
    U = []; LAMBDA = [];
    error('The matrix is not normal.');
  end

endfunction

//******************************************************************
// nulbasis  Basis for nullspace.
//
// N = nulbasis(A) returns a basis for the nullspace of A
// in the columns of N. The basis contains the n-r special 
// solutions to Ax=0.  freecol is the list of free columns.
//
// Example:
//
// >> A = [1 2 0 3;
//        [0 0 1 4];
//
// >> N = nulbasis(A)
//
//    N = [-2  -3]   
//        [ 1   0]
//        [ 0  -4]
//        [ 0   1]
//
// See also fourbase.
//******************************************************************
function N = nulbasis(A)

  // Surprisingly, Scilab accepts this statement (two results).
  [R, pivcol] = rref(A, sqrt(%eps));

  [m, n] = size(A);
  r = length(pivcol);
  freecol = 1:n;
  freecol(pivcol) = [];
  N = zeros(n, n-r);
  N(freecol, : ) = eye(n-r,n-r);
  N(pivcol,  : ) = -R(1:r, freecol);

endfunction

//******************************************************************
// orthcomp  Orthogonal complement of a subspace.
//
// BCOMP = orthcomp(B) returns a basis for the 
// orthogonal complement of the column space of B.
// This subspace contains all vectors orthogonal
// to the column space of B.
// It is the left nullspace of B.
//
// See also leftnull, nulbasis.
//******************************************************************
function BCOMP = orthcomp(B)

  BCOMP = leftnull(B);

endfunction

//******************************************************************
// partic  Particular solution of Ax=b.
//
// x = partic(A, b) returns a particular solution to Ax=b.
// This particular solution has all free variables set to zero.
// An empty vector is returned if Ax=b is not solvable.
//
// See also slash as in A\b .
//******************************************************************
function x = partic(A, b)

  [m, n] = size(A);
  [Rd, pivcol] = rref([A b]);
  r = length(pivcol);
  //
  // If the last column of the augmented matrix [A b] 
  // is a pivot column, then Ax=b has no solution.
  //
  if max(pivcol) == n+1
    x = [];
  else
    //
    // The values of the pivot variables are in the
    // last column (which is called d) of Rd.
    // The free variables are zero in this particular solution.
    //
    x = zeros(n, 1);
    d = Rd(:, n+1);
    x(pivcol) = d(1:r);
  end

endfunction

//******************************************************************
// plu  Rectangular PA=LU factorization *with row exchanges*.
//
// [P, L, U] = plu(A), for a rectangular matrix A, uses Gaussian elimination
// to compute a permutation matrix P, a lower triangular matrix L and 
// an upper trapezoidal matrix U so that PA = LU.
// U is the same size as A.  
// P and L are square, with as many rows as A.
// sign = det(P); it is 1 or -1.
//
// See also elim, slu, lu, rref, partic, nulbasis, colbasis.
//******************************************************************
function [P, L, U, pivcol, sign] = plu(A)

  [m, n] = size(A);
  P = eye(m, m);
  L = eye(m, m);
  U = zeros(m, n);
  pivcol = [];
  tol = sqrt(%eps);
  sign = 1;

  p = 1;
  for k = 1:min(m, n)
    [r, p] = findpiv(A(k:m, p:n), k, p, tol);
    if r ~= k
      A([r k], 1:n) = A([k r], 1:n);
      if k > 1
        L([r k], 1:k-1) = L([k r], 1:k-1);
      end // if
      P([r k], 1:m) = P([k r], 1:m);
      sign = -sign;
    end // if
    if abs(A(k, p)) >= tol
      pivcol = [pivcol p];
      for i = k+1:m
        L(i, k) = A(i, p) / A(k, p);
        for j = k+1:n
          A(i,j) = A(i, j) - L(i, k)*A(k, j);
        end // for
      end// for
    end // if
    for j = k:n
      U(k, j) = A(k, j) * (abs(A(k, j)) >= tol);
    end // for
    if p < n
      p = p+1;
    end // if
  end // for

  if argn(1) < 4
    nopiv = 1:n;
    nopiv(pivcol) = [];
    if ~isempty(pivcol)
      disp('Pivots in columns:')
      disp(pivcol);
    end
    if ~isempty(nopiv)
      disp('No pivots in columns:')
      disp(nopiv);
    end
    rank = length(pivcol);
    if rank > 0
      roworder = P*(1:m)';
      disp('Pivots in rows:')
      disp(roworder(1:rank)');
    end
  end

endfunction

//******************************************************************
// poly2str  Convert a polynomial coefficient vector to a string.
//
// p = poly2str(c) generates a string representation of the polynomial
// whose coefficents are in the vector c.  
// The default variable is 'x', unless otherwise specified by 
// p = poly2str(c, 's').
// The coefficients are approximated, if necessary, by the rational
// values obtained from rat.
//	
// If x has a numeric value and the elements of c are reproduced
// exactly by rat, then eval(poly2str(c)) will return the same value 
// as polyval(c, x).
//
// See also polyval, rat.
//******************************************************************
function p = poly2str(c)

  a = poly(c,"x","coef");
  p = pol2str(a);

endfunction

//******************************************************************
// project  Project a vector b onto the column space of A.
//
// p = project(A, b) returns the orthogonal projection of a 
// vector b onto the column space of A.
//
// [p, e] = project(A, b) also returns the vectors e = b - p.
// p is the projection of b onto the column space of A.
// e is the projection of b onto the left nullspace of A.
// Notice that b = p + e and p' * e = 0. 
//
// See also projmat.
//******************************************************************
function [p, e] = project(A, b)

  P = projmat(A);
  p = P * b;
  e = b - p;

endfunction
//******************************************************************
// projmat  Projection matrix for the column space of A.
//
// P = projmat(A) returns the projection matrix for
// the column space of A.
//******************************************************************
function P = projmat(A)

  A = colbasis(A);
  P = A*inv(A'*A)*A';

endfunction

//******************************************************************
// randperm  Random permutation.
//
// p = randperm(n) returns a random permutation of 1:n.
//******************************************************************
function p = randperm(n)

  p = fliplr(sort(rand(1, n)));

endfunction

//******************************************************************
// rowbasis  Basis for the row space. 
//
// B = rowbasis(A) returns a basis for the row space of A
// in the *columns* of B.
// The row space of A is the column space of A'.
// rowbasis finds the first r linearly independent 
// columns of A'.
//
// The command fourbase(A) uses rref(A) to find a 
// different basis for the row space of A.
//
// See also fourbase.
//******************************************************************
function B = rowbasis(A)

B   = colbasis(A');

endfunction

//******************************************************************
// samespan  Test if two matrices have the same column space.
//
// samespan(A1, A2) 
// If the column spaces of A1 and A2 are the same,
// the function returns the dimension of this subspace.
// If the subspaces are different, the function returns 0.
//
// See also rank.
//******************************************************************
function samespan(A1, A2)

  rankA1 = rank(A1)
  rankA2 = rank(A2)
  rankboth = rank([A1 A2])

  if (rankA1 == rankA2) & (rankA1 == rankboth) 
    disp('A1 and A2 have the same column space.');
  else
    disp('A1 and A2 have different column spaces.');
  end

endfunction

//******************************************************************
// signperm  Determinant of the permutation matrix with rows ordered by p.
//
// sign = signperm(p) returns the sign of the 
// permutation associated with the vector p.
//
// [sign, PERM] also returns the permutation matrix PERM.
//
// Example: Let p = [2 3 1].
// Then sign = 1 and PERM = [0 1 0; 0 0 1; 1 0 0] .
// 
//******************************************************************
function [sign, PERM] = signperm(p)

  n = length(p);
  I = eye(n,n);
  PERM = I(p, :);
  sign = det(PERM);

endfunction

//******************************************************************
// slu  LU factorization of a square matrix using *no row exchanges*.
//
// [L, U] = slu(A) uses Gaussian elimination to compute a unit 
// lower triangular L and an upper triangular U so that L*U = A.
// The algorithm will stop if a pivot entry is very small.
//
// See also slv, plu, lu.
//******************************************************************
function [L, U] = slu(A)

  [n, n] = size(A);

  for k = 1:n
    if abs(A(k, k)) < sqrt(%eps)
      disp(['Small pivot encountered in column ' int2str(k) '.'])
    end
    L(k, k) = 1;
    for i = k+1:n
      L(i,k) = A(i, k) / A(k, k);
      for j = k+1:n
        A(i, j) = A(i, j) - L(i, k)*A(k, j);
      end
    end
    for j = k:n
      U(k, j) = A(k, j);
    end
  end

endfunction

//******************************************************************
// x = slv(A, b)
//
// x = slv(A, b) tries to use the LU factorization computed 
// by slu(A) to solve the linear equation A*x = b.
// Since slu does *no row exchanges*, slv may fail if a 
// small pivot is encountered.
//
// See also slu, solve, slash, partic.
//******************************************************************
function x = slv(A, b)

  [L, U] = slu(A);

  // Forward elimination to solve L*c = b.
  // L is lower triangular with 1's on the diagonal.

  [n, n] = size(A);
  c = zeros(n, 1);
  for k = 1:n
  s = 0;
    for j = 1:k-1
      s = s + L(k, j)*c(j);
    end
    c(k) = b(k) - s;
  end

  // Back substitution to solve U*x = c.
  // U is upper triangular with nonzeros on the diagonal.

  x = zeros(n, 1);
  for k = n:-1:1
    t = 0;
    for j = k+1:n
      t = t + U(k, j)*x(j);
    end
    x(k) = (c(k) - t) / U(k, k);
  end
  x = x';

endfunction

//******************************************************************
// splu  Square PA=LU factorization *with row exchanges*.
//
// [P, L, U] = splu(A), for a square, invertible matrix A,
// uses Gaussian elimination to compute a permutation
// matrix P, a lower triangular matrix L and 
// an upper triangular matrix U so that P*A = L*U.
// P, L and U are the same size as A.
// sign = det(P); it is 1 or -1.
//
// See also slu, lu, rref, partic, nulbasis, colbasis.
//******************************************************************
function [P, L, U, sign] = splu(A)

  [m, n] = size(A);
  if m ~= n
    error('Matrix must be square.')
  end
  P = eye(n, n);
  L = eye(n, n);
  U = zeros(n, n);
  tol = sqrt(%eps);
  sign = 1;

  for k = 1:n
    if abs(A(k, k)) < tol
      for r = k:n
        if abs(A(r, k)) >= tol
          break
        end
        if r == n
          if argn(1) == 4
            sign = 0;
            return
          else
            disp('A is singular within tolerance.')
            error(['No pivot in column ' int2str(k) '.'])
          end
        end
      end
      A([r k], 1:n) = A([k r], 1:n);
      if k > 1, L([r k], 1:k-1) = L([k r], 1:k-1); end
      P([r k], 1:n) = P([k r], 1:n);
      sign = -sign;
    end
    for i = k+1:n
      L(i, k) = A(i, k) / A(k, k);
      for j = k+1:n
         A(i, j) = A(i, j) - L(i, k)*A(k, j);
      end
    end
    for j = k:n
      U(k, j) = A(k, j) * (abs(A(k, j)) >= tol);
    end
  end

  if argn(1) < 4
    roworder = P*(1:n)';
    disp('Pivots in rows:')
    disp(roworder');
  end

endfunction

//******************************************************************
// splv  The solution to a square, invertible system.
// x = splv(A, b) uses the PA = LU factorization
// computed by splu to solve Ax = b.
//
// See also slv, splu, slash.
//******************************************************************
function x = splv(A, b)

  [P, L, U] = splu(A);
  [n, n] = size(A);

  // Permute the right hand side.
  b = P*b;
         
  // Forward elimination to solve L*c = b.
  c = zeros(n, 1);
  for k = 1:n
    s = 0;
    for j = 1:k-1
      s = s + L(k, j)*c(j);
    end
    c(k) = b(k) - s;
  end

  // Back substitution to solve U*x = c.
  x = zeros(n, 1);
  for k = n:-1:1
    t = 0;
    for j = k+1:n
      t = t + U(k, j)*x(j);
    end
   x(k) = (c(k) - t) / U(k, k);
  end

endfunction

//******************************************************************
// symmeig  Eigenvalues and eigenvectors of a symmetric matrix.
// The matrix A is assumed to be symmetric.
//
// Q = symmeig(A) returns a set of orthonormal eigenvectors for A.
//
// [Q, LAMBDA] = symmeig(A) also returns the corresponding eigenvalues 
// on the diagonal of LAMBDA. The eigenvalues in LAMBDA are in 
// decreasing order.
//
// See also eigval, eigvec, eig.
//******************************************************************
function [Q, LAMBDA] = symmeig(A)

  [m, n] = size(A);
  if norm(A'-A) <= sqrt(%eps)
    //
    // The eigenvectors in S are linearly independent but not orthogonal.
    // Eigenvectors for different eigenvalues *are* orthogonal.
    // Gram-Schmidt (qr(S)) gives orthonormal eigenvectors in Q.
    //
    [S, LAMBDA] = spec(A);
    [Q, R] = qr(S);
  else
    Q = []; LAMBDA = [];
    error('The matrix is not symmetric.');
  end

endfunction

//******************************************************************
// tridiag  Tridiagonal matrix.
// T = tridiag(a, b, c, n) returns an n by n matrix that has 
// a, b, c as the subdiagonal, main diagonal, and superdiagonal 
// entries in the matrix.
//******************************************************************
function T = tridiag(a, b, c, n)

  T = b*diag(ones(n,1)) + c*diag(ones(n-1,1),1) + a*diag(ones(n-1,1),-1);

endfunction

