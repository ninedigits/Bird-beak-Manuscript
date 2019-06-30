function [angle, length] = birdbeak_length_angle2(g1, g2, a1, a2)
%{
Calculate the angle and length of endograft mal-apposition (ie birdbeaking)
using four coordinates (g1, g2, a1, a2), birdbeak height (bbh), and a
threshold that defines when the graft-aorta are close enough to be
considered touching. 

%}

angle = acos(dot((g1-g2),(a1-a2))/(norm(g1-g2)*norm(a1-a2)));
c = norm(g2 - a2);
length = norm(g1-g2);
angle = angle * (180/pi);

end

