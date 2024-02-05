//**********************************************************************
// File Name: C8.1-27.sci
//**********************************************************************

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct transformation matrix.
A = [2 1; 1 2];

// Construct angle vector.
theta = [0:2*%pi/50:2*%pi];

// Construct circle.
circle = [cos(theta); sin(theta)];

// Transform circle to ellipse.
ellipse = A * circle;

// Form aspect ratio.
square(-4,-4,4,4);

plot(circle(1,:),circle(2,:),ellipse(1,:),ellipse(2,:));

