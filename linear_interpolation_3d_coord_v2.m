function[smoothed_line_delta] = linear_interpolation_3d_coord_v2(smoothed_line,delta)

%% Generate interpolated line at delta arc length intervals
%
% This function takes an interpolated line with k points and generates a
% new interpolated line with m points, such that each adjacent point is
% exactly delta distance apart. 
%
% Inputs: 
%   1) smoothed_line: [nx3] matrix of coordinates that describe a curved line in 3D
%   space, where each coordinate is sequential in space; ie coord1 is
%   adjacent in space to coord2, and so on. 
%   2) delta: integer that describes newly desired spacing of interpolated line 
% Outputs:
%   1) smoothed_line_delta: [nx4] matrix where col1 describes total arc
%   length from proximal start and col2-4 describe the newly defined
%   coordinates derived from the original smoothed_line.

smoothed_line_length = length(smoothed_line);
smoothed_line_with_distances = zeros(smoothed_line_length, 1);
smoothed_line_with_distances = [smoothed_line_with_distances, smoothed_line];
line_distance = 0;
for i1=2:smoothed_line_length
    i0 = i1-1;
    coord1 = smoothed_line(i1,:);
    coord0 = smoothed_line(i0,:);
    d = norm(coord1-coord0);
    line_distance = line_distance+d;
    smoothed_line_with_distances(i1,1) = line_distance;
end
smoothed_line = smoothed_line_with_distances;

n_interp = round(smoothed_line(end,1))/delta;
smoothed_line_delta = zeros(n_interp + 1,4); 

inter_index = zeros(n_interp + 1,2);
smoothed_line_delta(1,:) = smoothed_line(1,:);
smoothed_line_delta(n_interp+1,:) = smoothed_line(end,:);
for i=2:n_interp
    smoothed_line_delta(i,1) = (i-1)*delta;
    for j=2:length(smoothed_line)
    	if((smoothed_line_delta(i,1) < smoothed_line(j,1)) && (smoothed_line_delta(i,1) > smoothed_line(j-1,1)))
            inter_index(i,1) = j-1; inter_index(i,2) = j; %inter_index saves the upper and lower bound index of the interpolation
        end
        if((smoothed_line_delta(i,1) == smoothed_line(j,1)))
            inter_index(i,1:2) = j;
            smoothed_line_delta(i,2) = smoothed_line(j,2);
        end
    end
    if(inter_index(i,1) ~= inter_index(i,2))
        % interpolation of x, y, z coord
        for j=2:4
            D1 = smoothed_line(inter_index(i,1),j);
            D2 = smoothed_line(inter_index(i,2),j);
            sx = smoothed_line_delta(i,1);
            s1 = smoothed_line(inter_index(i,1),1);
            s2 = smoothed_line(inter_index(i,2),1);
            Dx = D1+(D2-D1)*(sx-s1)/(s2-s1);
            smoothed_line_delta(i,j) = Dx;
        end
    end

%% Remove any remaining zero rows    
idx_zeros = smoothed_line_delta(:,2:4) ~= [0 0 0];        % index all coordinates with 3 zeroes
rows_to_keep = idx_zeros(:,1);                            % keep first row. tells us 
smoothed_line_delta = smoothed_line_delta(rows_to_keep,:);
end



 
 
 
 