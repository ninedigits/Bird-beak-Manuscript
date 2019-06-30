function [angle, curve] = angulation_curvature_at_point(smoothed_line,ref_point,window)

distance = window/2;
index_lz = find_index_of_array(ref_point,smoothed_line);
coord_lz = smoothed_line(index_lz,:);
coord_prox = find_proximal_point(smoothed_line,coord_lz,distance);
coord_dist = find_distal_point(smoothed_line,coord_lz,distance);
angle = angle_bw_3_points(coord_lz,coord_prox,coord_dist);
radius = radius_from_3_points(coord_prox,coord_lz,coord_dist);
curve = radius^-1;

end

