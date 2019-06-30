function [angle] = angle_bw_3_points(p1_center, p2, p3)

 
% Law of Cosines to determine angle between 3 points, with
% p1 as center point; ie the two lines intersect at p1


p12 = norm(p1_center-p2);
p13 = norm(p1_center-p3);
p23 = norm(p2-p3);

% Law of Cosines
cos_theta = (p12^2 + p13^2 - p23^2)/(2*p12*p13);
angle = (180/pi) * acos(cos_theta);
end








