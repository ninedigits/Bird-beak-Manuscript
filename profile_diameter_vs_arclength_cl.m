function [d_n, smoothed_markers, smoothed_arclength, interpolated] = ...
    profile_diameter_vs_arclength_cl(smoothed_cl_delta, aortic_contours)
%{

This function serves to show how the diameter of the aorta changes across
the center line of of the aorta.

    INPUTS:
        smoothed_cl_delta: [mx4] array where col1 is arclength and col2-4
        are x,y,z coordinates. eg: [0 x1 y1 z1; 0.1 x2 y2 z2; 0.2 x3 y3 z3 ...].
        the last coordinate
        
        aortic_contours: {1xn} cell array where each cell represents a
        contour S_i. 
    OUTPUTS:
        diameters_of_contour_i: actual diameter measurements; 1 for each
        contour. results in a [nx1].
        smoothed_markers: shows the markers on the smoothed line that
        corresponds to the closest location to the contours
        interpolated: 
    
    
    Algorithm Overview: 
    
    
    Given a centerline and n contours, we want to determine how diameter
    changes as a function of arclength from the proximal aorta to distal
    aorta.
    
    1) To begin, we need to know where centerline crosses through each
    contour. We estimate this finding the closest point on the centerline 
    to each respective contour centroid. Then measure the diameter of each
    contour, noting its location along the centerline. 
    
    2) Now that each diameter has a respective centerline marker, 
    we can count the distance along the centerline between each diameter to
    measure arclength. We use the arclength to interpolate the diameters. 
%}

smoothed_cl = smoothed_cl_delta(:,2:4); % xyz coordinates
arclength = smoothed_cl_delta(:,1);     % arclength
[~,n] = size(aortic_contours);          % n number of contours
smoothed_markers = zeros(n,3);          % mark the point on the cl
smoothed_indices = zeros(n,1);          % obtain index of marked point on cl
smoothed_arclength = zeros(n,1);        % store arclength
d_n = zeros(n,1);                       % store diameters of contours

%% 1) Determine diameters of contours and keep track along arc length
for i=1:n
    
    % Centroid of aortic contour. 
    centroid_Si = mean(aortic_contours{i}); 
    % Mark smoothed line and update variable
    smoothed_marker_i = find_closest_point(smoothed_cl, centroid_Si);
    smoothed_markers(i,:) = smoothed_marker_i;
    
    % Store indices of each marker and keep track of arclength
    index_i = find_index_of_row_in_array(smoothed_marker_i,smoothed_cl);
    smoothed_indices(i,:) = index_i;
    smoothed_arclength(i,:) = arclength(index_i,:);
    
    % Determine diameter of each contour for n number of contours
    d_n(i,:) = get_circumference(aortic_contours,i)/pi;
    
end 

%% interpolate diameters across arc length at 0.1 mm increments
interpolated = [];

syms d(x)

%diameter as a function of arclength to define interpolated values
count = 0;
for i=2:length(d_n)
    d2 = d_n(i,1);
    d1 = d_n(i-1,1);
    x2 = smoothed_arclength(i,1);
    x1 = smoothed_arclength(i-1,1);
    arcL = x1:0.1:x2; arcL = double(arcL)';
    L = round(x2-x1,1);
    x_in = 0:0.1:L; x_in = x_in';            
    d(x) = d1 + (x./L) * (d2 - d1);
    values_interp = double(d(x_in));
    interpolated = [interpolated; arcL, values_interp];
    
end



