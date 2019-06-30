clear all
close all
set(0,'DefaulttextInterpreter','none') %for plot3

%% Graft Sizing Data
graft_sizing = [40 34 45 31 34 34 28 34 34 31 37 34 37 40 45 45 37 31 40 45 31 40 40];

%%% Begin loop
for p=1:23
    
    %% Patient Data
    patients = {'T1', 'T2', 'T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20','T21','T22','T23'}; %all patients
    patient_ID = patients{p}
    all_modes = {'post_ID'};
    %mode_i = 'post_ID';
    
    %% Conditions for dissection vs aneurysm
    % Depending on whether or not the patient has a dissection or aneurysm,
    % there may not be aorta_fluid2.
    
    %file data: edit this when necessary for the appropriate file path
    file_path = "/Users/maxfrohlich/Dropbox/StanfordMatlabDevelopment/%s/"; 
    file_path = compose(file_path, patient_ID)
    
    
    %% Patient modes 
    for k=1:1
        
        aorta_name = 'aorta2';
        
        mode_i = all_modes{k};
        if k ~= 1 && p~=4
            continue
        end
        %loads aorta path data
        

        %% Load aorta contours
        contours_file = strcat('seg_', patient_ID, mode_i,'_',aorta_name) 
        contours = extract_contours(file_path, contours_file);
        scontours = struct;     % Initialize the structured array
        
        
        %% Pre-process contour files for graphing and calculations
        % Create structure file of contours for greenwich line code
        combined_contours = [];
        for i=5:length(contours)
            combined_contours = [combined_contours; contours{i}];
            field = 'coords';
            scontours(i).field = contours{i};
        end
        
        %% Determine centroids
        centroids = [];
        for i=1:length(contours)
            centroids = [centroids;mean(contours{i})];
        end
        
        %% Aorta Center of Mass
        
        % We use the center of mass to define the first coordinate of the 
        % inner and outer lines. Then we use the greenwich function to
        % perform vector projects across the contours.
        
        %%% Defining the aorta center of mass
        
        % Load the centerline. Then crop the centerline so that only the
        % arch is present. We define this by iterating through the
        % centerline until the reaching a z-coordinate value that is less
        % than the starting value. This tells us where the arch ends and 
        % where it begins to descend distally into the abdominal aorta. 
        % The purpose of this is to maintain consistency between the
        % different patients, since some models include more contours in
        % the abdominal aorta, which would otherwise 'weight' the center of 
        % mass to a more distal location.
        
        %cl_fn = strcat('coord_smooth_',patient_ID,mode_i,'_',aorta_name,'.txt');
        %cl = load(cl_fn);
        
        cl_z_crop = aorta_z_crop(centroids);
        Mzc = mean(cl_z_crop);      % Mean aorta centroid
        
        
        %% Find initial vector for the first aorta contour
        
        % Use this vector to define the first greenwich line coordinate,
        % which we use to determine the inner and outer lines of the aorta.
        
        temp_dist_array = [];
        s1 = contours{5};           %Skip first few aortic contours, since
                                    %some of them include the aortic valve
        for i=1:length(s1)
            % Record the distances of the first contour point to the
            % centroid ie center of mass
            temp_coord = s1(i,:);
            temp_dist = norm(temp_coord - Mzc(:,:));
            temp_dist_array = [temp_dist_array; temp_dist, temp_coord];
        end

        [min_inner, index] = min(temp_dist_array(:,1));
        [max_outer, indexo] = max(temp_dist_array(:,1));
        inner = temp_dist_array(index,2:4);
        outer = temp_dist_array(indexo,2:4);
        
        %% Plot how the first inner point is determined
        % This is solely for visualization purposes. It shows
        % how the first coordinate for the Greenwich line is 
        % determined. 
        
        inner_point_def = norm(inner-Mzc);
        
        t = inner_point_def;                    % Birdbeak height.
        t = 0:1:t;                              % Points for interpolation. 
        v = Mzc - inner;                        % Birdbeak vector.
        v = v/norm(v);                          % Unit vector of v.
        define_line = [inner;Mzc];
        
        for n = 1:length(t)-1
            define_line = [define_line; inner + n * v ];
        end
        
        
        %find the closest and furthest point from mean aorta centroid to establish
        %initial points


        
        
        %% Determine inner and outer lines using Greenwich line function

        % Initialize parameters
        innerpd = inner; 
        outerpd = outer;
        temp_iter = 5;
        length_Max = length(contours);
        gwline_inner = [];
        gwline_outer = [];
        %start at fiducial marker (ie the geometric center) and move from proximal to distal
        while temp_iter<length_Max-1
            %inner
            innerpd = find_green_point_for_next_contour(temp_iter,1,centroids,scontours,innerpd);
            gwline_inner = [gwline_inner; innerpd];
            
            %outer
            outerpd = find_green_point_for_next_contour(temp_iter,1,centroids,scontours,outerpd);
            gwline_outer = [gwline_outer;outerpd];
            
            %keep track of iterations
            temp_iter = temp_iter + 1;
        end
        
        % Combine files for easy plotting, however,
        % files are saved to disk independently for analysis
        gwline_aorta = [gwline_inner; gwline_outer]; % combined files
        
        %% Plot section
        
        figure
        plot3(gwline_inner(:,1),gwline_inner(:,2),gwline_inner(:,3),'.-','MarkerSize', 4,'LineWidth',3)
        hold on
        plot3(combined_contours(:,1),combined_contours(:,2),combined_contours(:,3),'.','Color',[0.49 0.18 0.56]),hold on
        %plot3(cgraft(:,1),cgraft(:,2),cgraft(:,3),'.','Color',[0.93 0.69 0.13]), hold on
        plot3(Mzc(1,1),Mzc(1,2),Mzc(1,3),'*'), hold on
        plot3(inner(1,1),inner(1,2),inner(1,3),'*'), hold on
        plot3(gwline_outer(:,1),gwline_outer(:,2),gwline_outer(:,3),'.-','MarkerSize', 4,'LineWidth',3),hold on
        plot3(define_line(:,1),define_line(:,2),define_line(:,3),'.-')
        ptitle = strcat('Patient ID: ', patient_ID);
        mtitle = strcat('Mode: ', mode_i);
        vtitle = strcat('Vessel name: ', aorta_name);
        title({ptitle, mtitle, vtitle});
        hold off
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis 'equal'
        
        %% Write files to disk
        cd('/Users/maxfrohlich/Dropbox/StanfordMatlabDevelopment/mainCodeFolder/inner_outer_data')
        fname_inner = strcat(patient_ID,mode_i,'_gwline_inner_',aorta_name,'.txt');
        fname_outer = strcat(patient_ID,mode_i,'_gwline_outer_',aorta_name,'.txt');
        save(fname_inner, 'gwline_inner', '-ascii')
        save(fname_outer, 'gwline_inner', '-ascii')
        cd('..')
        %% Graft inner and outer line. 
        % Muted for now because we are currently using a 
        % different method extracting the inner and outer lines.
        
        
        
%         %import graft data
%         if k==1
%             continue
%         else
%         graft_fn = strcat('seg_',patient_ID,mode_i,'_endograft2');
%         graft = extract_contours(file_path,graft_fn);
%         cgraft = [];
%         sgraft = struct;
%         gcfn = strcat(patient_ID,mode_i,'_endograft2_centroid.txt');
%         gcentroids = load(gcfn);
%         temp_iter = 1;
% 
% 
%         length_Max = length(graft);
%         g1 = graft{1};
%         for n=1:length(graft)
%             cgraft = [cgraft;graft{n}];
%             field = 'coords';
%             sgraft(n).field = graft{n};
%         end
%         %find closest and furthest distance to inital graft contour
%         temp_dist_array = []; %reset this variable for graft
%         for i=1:length(g1)-1
%             temp_coord = g1(i,:);
%             temp_dist = norm(temp_coord - Mzc(:,:));
%             temp_dist_array = [temp_dist_array; temp_dist, temp_coord];
%         end
%         [min_inner, index] = min(temp_dist_array(:,1));
%         [max_outer, indexo] = max(temp_dist_array(:,1));
%         ginner = temp_dist_array(index,2:4);
%         gouter = temp_dist_array(indexo,2:4);
%         ginnerpd = ginner;
%         gouterpd = gouter;
%         gwline_ginner = [ginner];
%         gwline_gouter = [gouter];
%         while temp_iter<length_Max-1
%             %inner
%             ginnerpd = find_green_point_for_next_contour(temp_iter,1,gcentroids,sgraft,ginnerpd);
%             gwline_ginner = [gwline_ginner; ginnerpd];
%             %outer
%             gouterpd = find_green_point_for_next_contour(temp_iter,1,gcentroids,sgraft,gouterpd);
%             gwline_gouter = [gwline_gouter;gouterpd];
% 
%             temp_iter = temp_iter + 1;
%         end
% 
%         %gwline_graft = [gwline_ginner; gwline_gouter];
% 
%         figure
%         plot3(gwline_inner(:,1),gwline_inner(:,2),gwline_inner(:,3),'bo')
%         hold on
%         plot3(gwline_outer(:,1),gwline_outer(:,2),gwline_outer(:,3),'ko')
%         hold on
%         plot3(combined_contours(:,1),combined_contours(:,2),combined_contours(:,3), '.')
%         hold on 
%         plot3(Mzc(:,1),Mzc(:,2),Mzc(:,3))
%         hold on
%         plot3(cgraft(:,1),cgraft(:,2),cgraft(:,3),'.')
%         hold off
%         ptitle = strcat('Patient ID: ', patient_ID);
%         mtitle = strcat('Mode: ', mode_i);
%         vtitle = strcat('Vessel name: ', aorta_name);
%         title({ptitle, mtitle, vtitle});
%         xlabel('x')
%         ylabel('y')
%         zlabel('z')
%         axis 'equal'
%         cd('/Users/maxfrohlich/Desktop/stanford_code/MATLAB_development_2017/Codes_original_max_edits/inner_outer_data')
%         fname_inner = strcat(patient_ID,mode_i,'_gwline_inner_',aorta_name,'.txt');
%         fname_outer = strcat(patient_ID,mode_i,'_gwline_outer_',aorta_name,'.txt');
%         fname_ginner = strcat(patient_ID,mode_i,'_gwline_ginner_endograft.txt');
%         fname_gouter = strcat(patient_ID,mode_i,'_gwline_gouter_endograft.txt');
%         dlmwrite(fname_inner,gwline_inner,' ')
%         dlmwrite(fname_outer,gwline_outer,' ')
%         %dlmwrite(fname_ginner,gwline_ginner,' ')
%         %dlmwrite(fname_gouter,gwline_gouter,' ')
%         cd('..')
%         end
    end
        
end


