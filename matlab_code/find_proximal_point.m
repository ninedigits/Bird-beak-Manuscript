function point = find_proximal_point(interp_line,starting_point,distance)

%% Find proximal point

% Take index of target point
% Then move backwards, toward the proximal end (ie toward the aortic root)

%% Starting index here of the point we want to analyze
index_of_start = find(ismember(interp_line,starting_point,'rows'));
dist = 0;                               % Keep track of distance
i = 0;                                  % Keep track of iteration
index_0 = index_of_start;       
index_1 = index_of_start - 1;           % Substract since we move backwards (ie dist -> prox)

%% Loop until target point found
try
    while dist<distance
        % Get two points and determine distance
        a = interp_line(index_0-i,:);       % Get first point
        b = interp_line(index_1-i,:);       % Get second point
        d = norm(a-b);                      % Get distance
        dist = dist + d;                    % Update distance
        i = i + 1;                          % Update iteration
    end
catch
end

%% Return target point
try
    point = a;                              % Return target point
catch
    point = starting_point;
end
end