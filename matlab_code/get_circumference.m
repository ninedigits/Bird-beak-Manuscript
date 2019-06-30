function C = get_circumference(contours, contour_num)

%{
Function takes a cell array of contours and determines the circumfrence 
of the contour given by contour_num
%}

C = 0;                                          % keep track of distance
my_contour = contours{contour_num};             % get specific contour
[total_iters, ~] = size(my_contour);            % find number of points in contour
% Duplicate the first point on the contour and add it to the end of the
% contour so that the final interpolation can be obtained.
my_contour = [my_contour; my_contour(1,:)];     
for i=1:total_iters
    p1 = my_contour(i,:);
    p2 = my_contour(i+1,:);
    d = norm(p1-p2);
    C = C + d; 
end

