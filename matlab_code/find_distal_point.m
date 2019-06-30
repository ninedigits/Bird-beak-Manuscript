function point = find_distal_point(interp_line, starting_point, distance_in_mm)
%% Find a coordinate along an interpolated line X distance from a starting coordinate
%
%
%
%%
index_of_start = find(ismember(interp_line,starting_point,'rows'));
index_of_start = index_of_start(1,1);

dist = 0;
i = 0;
index_0 = index_of_start;
index_1 = index_of_start + 1;
while dist<distance_in_mm
    a = interp_line(index_0+i,:);   % point 1
    b = interp_line(index_1+i,:);   % point 2
    d = norm(a-b);                  % distance b/w points
    dist = dist + d;                % update total distance
    i = i + 1;                      % update iter
end

point = a;
end