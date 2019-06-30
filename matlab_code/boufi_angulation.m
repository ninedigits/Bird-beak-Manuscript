function boufi_angle = boufi_angulation(smoothed_line,ref_point,window)

distance = window/2;
index_lz = find_index_of_array(ref_point,smoothed_line);
coord_lz = smoothed_line(index_lz,:);
coord_prox = find_proximal_point(smoothed_line,coord_lz,distance);
coord_dist = find_distal_point(smoothed_line,coord_lz,distance);

% goal is to get the instantaneous slop 
prox_ind = find_index_of_array(coord_prox,smoothed_line)-1; 
dist_ind = find_index_of_array(coord_dist,smoothed_line)+1; 

coord_prox_minus1 = smoothed_line(prox_ind,:); 
coord_dist_plus1 = smoothed_line(dist_ind,:); 

prox_vector = coord_prox_minus1 - coord_prox;
dist_vector = coord_dist_plus1 - coord_dist;

boufi_angle = angle_bw_2_planes(prox_vector, dist_vector);
end