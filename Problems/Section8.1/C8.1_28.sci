//**********************************************************************
// File Name: C8.1-28.sci
//**********************************************************************

//**********************************************************************
//
//  Name: plotArc
//
//  Purpose: The purpose of this function is to plot an arc.
//
//  Calling Sequence: plotArc(r,x,y,startAngle,endAngle)
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
//  Outputs:
//
//    None.
//
//**********************************************************************
function plotArc(r,x,y,startAngle,endAngle,angleResolution)

  // Contruct angle vector.
  theta = startAngle*%pi/180:2*%pi*angleResolution/180:endAngle*%pi/180;

  // Construct arc.
  arc = [cos(theta); sin(theta)];

  // Scale to the radius.
  arc = arc * r;

  plot2d(arc(1,:),arc(2,:));

endfunction


//**********************************************************************
//
//  Name: smile
//
//  Purpose: The purpose of this function is to plot a smile on a
//  circle.
//
//  Calling Sequence: smile(r)
//
//  Inputs:
//
//    r - The radius of curvature.
//
//    x - The x position of the smile.
//
//    y - The y position of the smile.

//  Outputs:
//
//    None.
//
//**********************************************************************
function smile(r,x,y)

  disp([r x y]);

  // Construct smile.
//  plotArc(r,x,y,260,280);

endfunction

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

plot(circle(1,:),circle(2,:));

// Draw smile.
//smile(0.7,0,0);


//plot(circle(1,:),circle(2,:),ellipse(1,:),ellipse(2,:));

