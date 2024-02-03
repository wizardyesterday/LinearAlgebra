//**********************************************************************
// File Name: C8.1-24.sci
//**********************************************************************

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct the house matrix.
H = [-6 -6 -7 0 7 6  6 -3 -3  0  0 -6;
     -7  2  1 8 1 2 -7 -7 -2 -2 -7 -7];

// Construct chimney matrix.
C = [5 5 4 4;
     3 6 6 4];

// Construct retrace matrix.
R = [ 6 6;
     -7 2];

// Construct the house plus chimney matrix.
HC = [H R C];

// Plot the complete matrix.
plot(HC(1,:),HC(2,:));


