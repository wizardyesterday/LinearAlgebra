//**********************************************************************
// File Name: C8.1-28.sci
//**********************************************************************

//**********************************************************************
//
//  Name: plotArc
//
//  Purpose: The purpose of this function is to plot an arc.  Notice
//  the order of operations:
//    1. Scale the radius of the arc.
//    2. Transform the arc by the transformation matrix.
//    3. Translate the arc to the appropriate position.
//    4. Draw the arc.
//
//  Calling Sequence: plotArc(r,x,y,startAngle,endAngle,A)
//
//  Inputs:
//
//    r - The radius of curvature.
//
//    x - The x reference position of the arc.
//
//    y - The y reference position of the arc.  Together (x,y) define the
//    origin.
//
//    startAngle - The start angle of the arc, in units of degrees.
//
//    endAngle - The end angle of the arc in units of degrees.
//
//    angleResolution - The resolution of the arc in units of degrees.
//
//    A - The transformation matrix.
//
//  Outputs:
//
//    None.
//
//**********************************************************************
function plotArc(r,x,y,startAngle,endAngle,angleResolution,A)

  // Contruct angle vector.
  theta = startAngle*%pi/180:2*%pi*angleResolution/180:endAngle*%pi/180;

  // Construct arc.
  arc = [cos(theta); sin(theta)];

  // Compute matrix form.
  s = size(arc);

  // Create a matrix of 1's with the same structure of the arc matrix.
  o = ones(s(1),s(2));

  // Scale for translation of the arc.
  o(1,:) = o(1,:) * x;
  o(2,:) = o(2,:) * y;

  // Scale to the radius.
  arc = arc * r;

  // Transform.
  arc = A * arc;

  // Translate.
  arc = arc + o;

  // Draw the arc.
  plot2d(arc(1,:),arc(2,:));

endfunction

//**********************************************************************
//
//  Name: plotSmiley
//
//  Purpose: The purpose of this function is to plot a smile on a
//  circle.
//
//  Calling Sequence: plotSmiley(r,x,y,A)
//
//  Inputs:
//
//    r - The radius of curvature.
//
//    x - The x reference position of the smile.
//
//    y - The y reference position of the smile.  Together (x,y) define the
//    origin.
//
//    A - The transformation matrix.
//  Outputs:
//
//    None.
//
//**********************************************************************
function plotSmiley(r,x,y,A)

  // Construct smile.
  plotArc(r,x,y,250,290,1,A);

endfunction

//**********************************************************************
//
//  Name: plotEyes
//
//  Purpose: The purpose of this function is to plot eyes on a
//  circle.
//
//  Calling Sequence: plotEyes(r,x,y,A)
//
//  Inputs:
//
//    r - The radius of curvature.
//
//    x - The x reference position of the smile.
//
//    y - The y reference position of the smile.  Together (x,y) define the
//    origin.
//
//    A - The transformation matrix.
//  Outputs:
//
//    None.
//
//**********************************************************************
function plotEyes(r,x,y,A)

  // Compute coordinates for left eye.
  xLeft = 0.3 * cos(5*%pi/6);
  yLeft = 0.3 * sin(5*%pi/6);

  // Compute coordinates for left eye.
  xRight = 0.3 * cos(%pi/6);
  yRight = 0.3 * sin(%pi/6);

  // Plot left eye.
  plotArc(r,xLeft,yLeft,0,360,0.1,A);

  // Plot right eye.
  plotArc(r,xRight,yRight,0,360,0.1,A);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct ellipse transformation matrix.
A = [2 1; 1 2];

// Construct identity transformation matrix.
I = eye(2,2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Draw smiley face.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(211);
title('Untransformed');
isoview(-1,1,-1,1);

// Construct circle, referenced to the origin.
plotArc(1,0,0,0,360,1,I);

// Draw smile, referenced to the origin.
plotSmiley(0.7,0,0,I);

// Draw eyes.
plotEyes(0.05,0,0,I);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Draw transformed smiley face.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(212);
title('Transformed With [2 1; 1 2]');
isoview(-3,3,-3,3);

// Construct circle, referenced to the origin.
plotArc(1,0,0,0,360,1,A);

// Draw smile, referenced to the origin.
plotSmiley(0.7,0,0,A);

// Draw eyes.
plotEyes(0.05,0,0,A);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


