//**********************************************************************
// File Name: C3.5_27.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct checkerboard rows.
row1 = [1 0 1 0 1 0 1 0]
row2 = [0 1 0 1 0 1 0 1];

// Construct checkeboard matrix.
B = [row1; row2; row1; row2; row1; row2; row1; row2];

// Compute the left nullspace basis.
B_left_null = nulbasis(B');

// Choose values for the chess pieces as used by the Sargon program.
pawn = 1;
knight = 2;
bishop = 3;
rook = 4;
queen = 5;
king = 6;

// Construct rows of the chess board.
outer = [rook knight bishop queen king bishop knight rook];
inner = pawn * ones(1,8);
zero = zeros(4,8);

// Construct matrix representation of chess board.
C = [outer; inner; zero; inner; outer];

// Compute the left nullspace basis.
C_left_null = nulbasis(C');

// Compute the nulspace basis.
C_null = nulbasis(C);

