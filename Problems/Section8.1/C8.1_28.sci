//**********************************************************************
// File Name: C8.1-28.sci
//**********************************************************************

//**********************************************************************
//
//  Name: plotArc
//
//  Purpose: The purpose of this function is to plot an arc.
//
//  Calling Sequence: plotArc(r,x,y,startAngle,endAngle,A)
//
//  Inputs:
//
//    r - The radius of curvature.
//
//    x - The x position of the arc.
//
//    y - The y position of the arc.
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

  plot2d(arc(1,:),arc(2,:));

endfunction

//**********************************************************************
//
//  Name: smile
//
//  Purpose: The purpose of this function is to plot a smile on a
//  circle.
//
//  Calling Sequence: smile(r,x,y,A)
//
//  Inputs:
//
//    r - The radius of curvature.
//
//    x - The x position of the smile.
//
//    y - The y position of the smile.
//
//    A - The transformation matrix.
//  Outputs:
//
//    None.
//
//**********************************************************************
function smile(r,x,y,A)

  disp([r x y]);

  // Construct smile.
  plotArc(r,x,y,250,290,1,A);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct ellipse transformation matrix.
A = [2 1; 1 2];

// Construct identity transformation matrix.
I = eye(2,2);

// Form aspect ratio.
square(-4,-4,4,4);

// Construct circle.
plotArc(1,1,1,0,360,1,I);

// Draw smile.
smile(0.7,1,1,I);
