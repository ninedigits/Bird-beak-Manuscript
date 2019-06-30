function [angles_all, curves_all, x_ang_curv] = curvature_and_angulation_with_arclength(smoothed_line_delta, window_min, window_max, window_step)
global ax
%{

This function takes an interpolated line and determines the curvature and
angulation varied by a different windows sizes. The goal of this function
will help us understand if there is an optimal curvature or angulation
window size to use.

PARAMETERS:

Input
    smoothed_line_delta: a 3d line described by an [nx4] matrix where col1
    is the current arclength, where:
        smoothed_line_delta(:,1) = 0
        smoothed_line_delta(:,2) = 0.1
        smoothed_line_delta(:,2) = 0.2
        ...

    window_min: the smallest size of the window, where each unit is a 1 mm
    window_max: the maximum size of the window
    window_Step: the step size of each window iteration
    Example:    window_min = 10
                window_max = 14
                window_step = 2
                outputs: [14 12 10] windows sizes.
Ouputs
    curve_angulation: [nx3] array where:
        col1: arclength distance from proximal end
        col2: curvature
        col3: arclength
%}

smoothed_line = smoothed_line_delta(:,2:4);
delta = smoothed_line_delta(:,1);
windows_all = zeros(1,(window_max-window_min)/window_step);
legend_label = {};

%% Establish windows
for i=0:length(windows_all)
    windows_all(1,i+1) = window_max - i*window_step;
    legend_label = [legend_label, num2str(window_max - i*window_step)];
end

%% find initial coordinate
p_mid = find_distal_point(smoothed_line,smoothed_line(1,:),window_max/2);
iter = find_index_of_array(p_mid, smoothed_line);
%%

step_size = 1;
angles = [];
curves = [];
coordinates = [];
curves_all = [];
angles_all = [];
try
    loop_on = 1;
    while loop_on
        %Need to start at midpoint of max window size, then move back a
        %window length that depends on which window size is currently
        %selected.
        angles_while_loop = [];
        curves_while_loop = [];
        
        % Starting midpoint with max windows size;
        p_mid = smoothed_line(iter,:);
        for i=1:length(windows_all) %% For each window size:
            w = windows_all(1,i)/2; %% We use half b/c mid to prox is 1/2 
            % the window
            
            % Prox and Dist points start at max window size and then
            % decrease. eg 14, 12, 10 correspond to the window sizes where
            % w = 1/2 * window_size.
            p_prox = find_proximal_point(smoothed_line,p_mid,w);
            p_dist = find_distal_point(smoothed_line,p_mid,w);
            
            angles_while_loop = [angles_while_loop, angle_bw_3_points(p_mid,p_prox,p_dist)]; %, p_prox, p_mid, p_dist]; 
            curves_while_loop = [curves_while_loop, radius_from_3_points(p_prox,p_mid,p_dist)^-1];
        end
        arclength_iter = delta(iter,1);
        angles_all = [angles_all; arclength_iter, angles_while_loop];
        curves_all = [curves_all; arclength_iter, curves_while_loop];
        
        % Update Iter
        iter = iter + step_size;
        
    end
catch
    loop_on = 0;

end

x_ang_curv = angles_all(:,1);
y_angles = angles_all(:,2:end);
y_curves = curves_all(:,2:end);






end