function marker = find_closest_point(interp_line, point)

%{
This function takes an interpolated line and a point, and then finds the
closest point on the interpolated line, called the marker, using the point.

INPUTS: 
    1) interp_line: [nx3] matrix that forms a line of points.
    2) point: [1x3] matrix that is a coordinate
OUTPUT:
    1) marker: the point on the interpolated line that is closest to the
    reference point.

%}
if ~isnumeric(interp_line)
    type = class(interp_line);
    error(compose("variable: 'interp_line' inputs must be a numberical [nx4] array; instead got %s",type))
elseif ~isnumeric(point)
    type = class(point);
    error(compose("variable: 'point' inputs must be a numberical [1x3] array; instead got %s", type))
end


[line_size, ~] = size(interp_line);
dist_coords = zeros(line_size,4); %col1 = distance col2-4=coordinate 
for i=1:line_size
    p = interp_line(i,:);
    d = norm(p-point(1,:));
    dist_coords(i,:) = [d,p];    
end

[~, ind] = min(dist_coords(:,1));
marker = dist_coords(ind,2:4);