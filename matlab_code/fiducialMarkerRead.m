function fiducial_output = fiducialMarkerRead(arclength, smoothed_line, LCCA_inner_curve_marker)
%{
This function translates an arclength along an interpolated line with a
given starting point.

Args:
    arclength (int): The arclength distance to travel. Can be positive or
    negative. Positive values indicate a distal direction, negative values
    indicate a proximal direction.
 
    smoothed_line (matrix[nx3]): a matrix of coordinates that represents an
    interpolated line. cols 1, 2, 3 --> x, y, z.

    LCCA_inner_curve_marker (matrix[1x3]): the point on the smoothed_line that is 
    closest to the LCCA ostium. cols 1, 2, 3 --> x, y, z.

Returns:
    fiducial_output (matrix[1x3]): the resulting point after translating
    the arclength distance away from LCCA_marker.
%}

if arclength > 0
    fiducial_output = find_distal_point(smoothed_line, LCCA_inner_curve_marker, arclength);
elseif arclength < 0
    
    fiducial_output = find_proximal_point(smoothed_line, LCCA_inner_curve_marker, abs(arclength));
else
    fiducial_output = LCCA_inner_curve_marker;
end


end 