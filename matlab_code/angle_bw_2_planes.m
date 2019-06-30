function angle = angle_bw_2_planes(n1, n2)
%{
Function takes normal vector of two planes, n1 and n2,
and returns the angle between these planes. 
%}

angle = (180/pi) * acos(abs(dot(n1,n2))/(norm(n1)*norm(n2)));

end