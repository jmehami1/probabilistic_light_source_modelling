function [ligLoc, ligDir] = SolveLightSrcDirLoc(ptReflFrame,normRefl, numPts)
%Solves the triangulation problem of a point source's location and
%principal direction. This light source is with a frame
%camera where the frame camera's coordinate frame is the world coordinates.

% INPUTS:
%       ptRelfFrame - brightest points found on a reflective sphere in
%       different poses. Matrix numPts x 3
%       normRefl - the surface normal at each brightest point on the
%       reflective sphere. Matrix of numPts x 3
%       numPts - number of brightest points. Number of different poses.
% OUTPUTS:
%       ligLoc - the point light source's location relative to the frame
%       camera
%       ligDir - the principal direction of the point light source

%Author: Jasprabhjit Mehami, 13446277

%A matrix
A = zeros(3*numPts, 3);
%b matrix
b = zeros(3*numPts, 1);

for i = 1:numPts
    %components of the point on the reflective sphere
    x1 = ptReflFrame(i,1);
    x2 = ptReflFrame(i,2);
    x3 = ptReflFrame(i,3);
    
    %components of the surface normal at the poi on reflective sphere
    n1 = normRefl(i,1);
    n2 = normRefl(i,2);
    n3 = normRefl(i,3);
    
    A(i, :) = [(1 - 2*n1^2), (-2*n1*n2), (-2*n1*n3)];
    A(numPts + i, :) = [(-2*n1*n2), (1 - 2*n2^2), (-2*n2*n3)];
    A(2*numPts + i, :) = [(-2*n1*n3), (-2*n2*n3), (1 - 2*n3^2)];
    
    b(i) = -(2*n1^2*x1 - 2*x1 + 2*n1*n2*x2 + 2*n1*n3*x3);
    b(numPts + i) = -(2*n2^2*x2 - 2*x2 + 2*n1*n2*x1 + 2*n2*n3*x3);
    b(2*numPts + i) = -(2*n3^2*x3 - 2*x3 + 2*n1*n3*x1 + 2*n2*n3*x2);
end

%solve the least squares problem to find the source location
[ligLoc, conFlag] = lsqr(A, b);

if conFlag
   disp('Light location did not converge');
end

for i = 1:numPts
    %components of the point on the reflective sphere
    x1 = ptReflFrame(i,1);
    x2 = ptReflFrame(i,2);
    x3 = ptReflFrame(i,3);
    
    %components of the surface normal at the poi on reflective sphere
    n1 = normRefl(i,1);
    n2 = normRefl(i,2);
    n3 = normRefl(i,3);
    
    s1 = ligLoc(1);
    s2 = ligLoc(2);
    s3 = ligLoc(3);
    
    q = sqrt(s1^2 - 2*s1*x1 + s2^2 - 2*s2*x2 + s3^2 - 2*s3*x3 + x1^2 + x2^2 + x3^2);
    
    A(i, :) = [(- q + n1^2*q*2), (n1*n2*q*2), (n1*n3*q*2)];
    A(numPts + i, :) = [(n1*n2*q*2), (- q + n2^2*q*2), (n2*n3*q*2)];
    A(2*numPts + i, :) = [(n1*n3*q*2), (n2*n3*q*2), (- q + n3^2*q*2)];
    
    b(i) = x1;
    b(numPts + i) = x2;
    b(2*numPts + i) = x3;
end

%solve the least squares problem to find the principal direction
[ligDir, conFlag] = lsqr(A, b);

if conFlag
   disp('Principal direction did not converge');
end

%normalise direction vector
ligDir = ligDir./norm(ligDir);

end
