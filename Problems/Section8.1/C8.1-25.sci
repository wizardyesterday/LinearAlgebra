//**********************************************************************
// File Name: C8.1-25.sci
//**********************************************************************

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct the house matrix.
H = [-6 -6 -7 0 7 6  6 -3 -3  0  0 -6;
     -7  2  1 8 1 2 -7 -7 -2 -2 -7 -7];

// Construct linear transformation matrices.
A1 = [1 0; 0 1];
A2 = [cos(35*%pi/180) -sin(35*%pi/180); sin(35*%pi/180) cos(35*%pi/180)];
A3 = [0 1; 1 0];
A4 = [0.7 0.3; 0.3 0.7];

// Perform transformations.
H11 = A1 * H;  H12 = A1' * A1 * H;
H21 = A2 * H;  H22 = A2' * A2 * H;
H31 = A3 * H;  H32 = A3' * A3 * H;
H41 = A4 * H;  H42 = A4' * A4 * H;

// Plot results.
subplot(241);
title("A1 * H");
plot(H11(1,:),H11(2,:));
subplot(242);
title("A1t * A1 * H");
plot(H12(1,:),H12(2,:));

subplot(243);
title("A2 * H");
plot(H21(1,:),H21(2,:));
subplot(244);
title("A2t * A2 * H");
plot(H22(1,:),H22(2,:));

subplot(245);
title("A3 * H");
plot(H31(1,:),H31(2,:));
subplot(246);
title("A3t * A3 * H");
plot(H32(1,:),H32(2,:));

subplot(247);
title("A4 * H");
plot(H41(1,:),H41(2,:));
subplot(248);
title("A4t * A4 * H");
plot(H42(1,:),H42(2,:));

