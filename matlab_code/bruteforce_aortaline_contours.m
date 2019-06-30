function graft_inner_points = bruteforce_aortaline_contours(aortaline,graft_contours,spacing)

%{
This function takes an interpolated line and a cell array, where each cell
is an contour of points. Using a bruteforce approach, this function
determines the minimum point within each cell (ie contour) that lies 
closest to the points on the interpolated line. 

INPUT
    aortaline:      [nx3] matrix that describes the inner line of the aorta.
    graft_contours: {1xm} cell array, where each m is an nx3 matrix of
                    contour points.
    spacing:        skip x amount of points when itering through the aorta
OUTPUT
    graft_inner_points: an [nx4] matrix that describes the closest point on
                        contour to the aorta inner line. Columns 1:3
                        describe coordinate location; column 4 is the
                        min distance between the line and the contour

%}
fprintf('\nextracting graft innerline...\n')
graft_inner_points = [];
n = length(graft_contours);
upd = textprogressbar(n);
for i=1:length(graft_contours)
    upd(i);
    contour_distances = [];
    %look through the cell array
    Gi = graft_contours{i}; %ith contour
    all_coords = [];
    for j=1:length(Gi) %for each point within a contour
        j;
         %for storing distances and coords
        %look through each point on contour
        Cj = Gi(j,:); %jth point of the the ith contour
        d = Inf;
        coords = [];
        for k=1:spacing:length(aortaline) %larger spacing reduces comp time
            k;
            Ak = aortaline(k,:); %kth point on the aorta innerline
            dk = norm(Ak-Cj);
            if dk < d
                d = dk; %update d when distance is shorter than all previous
                coords = Cj; %need to update coords
            end
            
            
        end
        %coords
        %d
        all_coords = [all_coords; coords, d];

        
    end
    [~, index] = min(all_coords(:,4));
    graft_point = all_coords(index,:);
    graft_inner_points = [graft_inner_points;graft_point]; 
    
end

end