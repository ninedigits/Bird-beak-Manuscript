function line = reference_line(p1,p2)
% takes two points and returns an interpolated line b/w the points


d = norm(p1-p2); % distance b/w points

                            
t = 0:1:d;                  % Points for interpolation. 
v = p2 - p1;                % Vector
v = v/norm(v);              % Unit vector of v.

[~,n] = size(t);
line = zeros(n,3);          % Initialize points
line(1,:) = p1;
line(end,:) = p2;
for n = 1:length(t)-1
    line(n,:) = p1 + n * v ;
end
end