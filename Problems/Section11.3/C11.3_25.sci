//**********************************************************************
// File Name: C11_3.25.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Functions.
//**********************************************************************
//
//**********************************************************************
//
//  Name: CG_IterationKludge
//
//  Purpose: The purpose of this function is to solve Ax = b using the 
//  conjugate iteration.
//
//  Calling Sequence: [x,x10,x20] = CG_IterationKludge(A.b)
//
//  Inputs:
//
//    A - The positive definite input matrix.
//
//    b - The input vector used to set initial conditions/
//
//  Outputs:
//
//    x - The solution to the equation Ax = b.
//
//**********************************************************************
function [x,x10,x20] = CG_IterationKludge(A,b)

  // Compute order of A.
  N = size(A);
  N = N(1);

  // Enforce column  vector.
  b = b(:);

  // Initial conditions.
  x0 = zeros(N,1);
  r0 = b;
  d0 = r0;
    
  for n = 1:N

    // Step length x_n-1 to x_n.
    Alpha = (r0' * r0) / (d0' * A * d0);

    // Approximate solution.
    x = x0 + Alpha * d0;

    //New residual b - A * x_n.
    r = r0 - Alpha * A * d0;

    if n == 10
      x10 = x;
    else
      if n == 20
        x20 = x;
      end
    end

    // Improve this step.
    Beta = (r' * r) / (r0' * r0);

    // Next search iteration.
    d = r + Beta * d0;

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // We need a way to bail out of this loop, otherwise,
    // when a solution is found, th next iteration will
    // produce a division by zero error.  The algorithm,
    // stated in the textbook exercises, really should
    // have included this bail out clause.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    if norm(r) == 0
      break
    else
      // Update x_n-1, r_n-1, d_n-1.
      x0 = x;
      r0 = r;
      d0 = d;
    end
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Create -1,2,-1 matrix.
A = secondDifferenceMatrix(100);;

// Create (1,...,1) columnvector.
b = ones(100,1);

//Compute solution to Ax=b.
[x,x10,x20] = CG_IterationKludge(A,b);

// Compute exact solution.
x_exact = A \ b;

// Plot results.
subplot(411);
title('Solution at 10thc Iteration');
plot(1:100,x10);

subplot(412);
title('Solution at 20thc Iteration');
plot(x20);

subplot(413);
title('Solution at Full Iteration');
plot(x);

subplot(414);
title('Exact solution');
plot(x_exact);



