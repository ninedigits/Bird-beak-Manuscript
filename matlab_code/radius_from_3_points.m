function R = radius_from_3_points(a,b,c)
%{
This function determines the radius of a circle defined by 3 points

INPUTS:
    a: first point
    b: second point
    c: third point

OUTPUTS:
    R: radius of curvature
%}

%test case
% a = [2 0 0];
% b = [0 2 0];
% c = [-2 0 0];

%length of lines formed
AB = norm(a-b);
AC = norm(a-c);
BC = norm(b-c);

%Heron's formula to determine area of triangle given 3 points
s = 1/2 * (AB + AC + BC);
K = sqrt(s*(s-AB)*(s-AC)*(s-BC)); %area of triangle

%Radius of curvature
R = AB * AC * BC / (4*K);