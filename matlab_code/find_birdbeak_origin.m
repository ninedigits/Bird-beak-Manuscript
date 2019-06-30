function op = find_birdbeak_origin(aorta_line,graft_line, threshold)
%{
This function finds the birdbeak origin and outputs the arclength 
and two origin points: the graft origin and corresponding aorta origin.

INPUT:
    aorta_line: the smoothed aorta line. nx3
    graft_line: the splined graft line. nx3
    threshold: the cutoff value when we say that the graft and aorta touch.
    

OUTPUT:
    op: 1x7 matrix in form [arclength, graft origin, aorta origin]. 'op'
    stands for origin points.

        arclength:op(1,1) the distance between the first graft inner point
        and the contact point on graft.
        
        graft origin: op(1,2:4) the point on the graft that travels an
        arclength from the first graft point until the origin point

        aorta origin: op(1,5:7) the corresponding aorta point that is
        closest to the graft origin.
        

%}

fprintf('\nFinding endograft-aorta contact location.......')
graft_line = graft_line(:,1:3);
length(aorta_line);

arclength = 0;

for i=1:length(graft_line) %for each point on graft (~120 points per graft)
    
    g = graft_line(i,:);
    d = Inf; %keep the smallest distance
    c = [0 0 0]; %update this point when smallest distance found
    %now cycle through every point on aorta and find the closest point
    for j=1:length(aorta_line)
        a = aorta_line(j,:);
        %calculate distance between points
        temp_d = norm(g - a);

        
        if temp_d < d
            d = temp_d; %update distance
            c = a;      %update coordinates on aorta
            
        elseif d < temp_d
            continue %continue to next iter if the distance is increase 
        end
        %stop iterating when the distance begins increasing to save
        %resources
    end
    if i>1
        arclength = arclength + norm(g-graft_line(i-1,:));
    end
   
        
    if d<=threshold && i>2 %needs to go at least one iter 
        break
    end
   
end

op = [arclength, g, c];
fprintf('\nCompleted |..................................................| Done.')
end

