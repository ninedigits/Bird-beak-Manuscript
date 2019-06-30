function fiducialArclength = fiducialMarkerWrite(LCCA_inner_curve_marker, smoothed_line, target)
%{
This function can translate a point on a post-operative model to a
pre-operative model. It takes a line and two points, and determines the
arclength distance to those points. Arclength can be positive 
(prox  --> distal) or negative (prox <-- distal).

Args:
    LCCA_inner_curve_marker (matrix[1x3]): the point on the smoothed_line
    where the LCCA ostium, is closest to 
%}
LCCA_innerline_marker_index = find_index_of_row_in_array(LCCA_inner_curve_marker, smoothed_line);
landing_zone_index = find_index_of_row_in_array(target, smoothed_line);

%%
if landing_zone_index >= LCCA_innerline_marker_index
    % then the LZ is distal to the LCCA marker
    orientation = 1;
    smoothed_line_segment = smoothed_line(LCCA_innerline_marker_index:landing_zone_index,:);
else landing_zone_index < LCCA_innerline_marker_index;
    % then the LZ is proximal to the LCCA marker
    orientation = -1;
    smoothed_line_segment = smoothed_line(landing_zone_index:LCCA_innerline_marker_index,:);
end

fiducialArclength = orientation * arclength_distance(smoothed_line_segment);
end