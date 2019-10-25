%% Birdbeak Analysis

%{

Author:         Max Frohlich
Date:           Nov 2017
Affliations:    Stanford University and San Jose State University.
E-mail:         max.frohlich@gmail.com

%}


clear all
close all
%%
%%
%% Username

username = 'maxfrohlich';
%username = 'gysuh';

%% Save data 
% Use 0 for testing
Save_num = 0; %16 current for ID testing plp2 metric 15 Feb 2018



tic     % timer
%% Patient data
patients = {'T1', 'T2', 'T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20','T21','T22','T23'}; %all patients
disorder_type = ['D' 'D' 'A' 'D' 'D' 'D' 'A' 'D' 'A' 'A' 'D' 'D' 'D' 'A' 'D' 'D' 'A' 'D' 'A' 'D' 'A' 'A' 'A']';
graft_sizing = [40 34 45 31 34 34 28 34 34 31 37 34 37 40 45 45 37 31 40 45 31 40 40];
all_modes = {'post_ID','pre_ID'};


%% Adjustable parameters:

threshold = 3;      %  Endograft-Aorta Contact Threshold

%% Patient IDs and output data
%patient_IDs = [1:7,9:23];
patient_IDs = 9;
[~, total_patients] = size(patient_IDs);
preop_geometric_output_data = zeros(total_patients,6);
postop_geometric_output_data = zeros(total_patients,6);
birdbeak_and_graft = zeros(total_patients,5);

graft_expansion_index = [];
boufi_angles = [];
bbh_all = [];

patient_count = 1;
for p=patient_IDs % cycle through all patients.
    
    
    patient_ID = patients{p};    
    file_path = "/Users/%s/Dropbox/StanfordMatlabDevelopment/%s/";
    file_path = compose(file_path,username,patient_ID);
    
    
    % Cycle through all modes
    output_LZ_preop = zeros(1,4);
    figure
    for k=1:length(all_modes)
                
        mode_i = all_modes{k}; %current mode
        
        % Patient 4 has overlapping contours in the 'aorta2' file
        if p==4
            aorta_name = 'aorta_fluid2';
        else
            aorta_name = 'aorta2';
        end

        fprintf(strcat('\nPatient: \t\t',patients{p},'\nMode: \t\t\t',mode_i,'\nAorta mode: \t\t',aorta_name,'\n'));
        fprintf('\n\n')
        
        %% Aorta Load Contours
        % Load aortic contours files
        % Keep original contour file as cell array 
        % Create additional matrix from cell array contents for plotting
        
        contours_fn = strcat('seg_', patient_ID, mode_i,'_',aorta_name); 
        contours_aorta = extract_contours(file_path, contours_fn);
        
        contours_aorta_s = length(contours_aorta);
        combined_contours_aorta = []; 
        centroids_aorta = [];
        % Loops will change size on each patient
        for i=1:contours_aorta_s            
            %combine contours into one matrix for plotting 
            combined_contours_aorta = [combined_contours_aorta;contours_aorta{i}];
            %combine centroids into one matrix for plotting 
            centroids_aorta = [centroids_aorta; mean(contours_aorta{i})];
        end
       
        %% Import Smoothed Aorta Innerline. 

        
        % 2 paths for different file locations
        pathname1 = compose("/Users/%s/Dropbox/StanfordMatlabDevelopment/smooth_curve_data/",username);
        pathname2 = compose("/Users/%s/Dropbox/StanfordMatlabDevelopment/%s",username,patient_ID);
        
        % Smoothed inner line of aorta
        smoothed_line_fp = compose("%scoord_smooth_%s%s_gwline_inner_%s.txt",pathname1,patient_ID,mode_i,aorta_name);        
        smoothed_line = load(smoothed_line_fp);
        
        % Smoothed centerline
        smoothed_cl = load(compose("%s/coord_smooth_%s%s_%s.txt",pathname2,patient_ID,mode_i,aorta_name));
        %% Downsample interpolation
        % Downsample smoothed_line in delta-distance increments
        
        delta = 0.1;
        smoothed_line_delta = linear_interpolation_3d_coord_v2(smoothed_line, delta);
        smoothed_line = smoothed_line_delta(:,2:4); %col1 is distances, which we dont need for plotting
        
        % Repeat for centerline
        smoothed_cl_delta = linear_interpolation_3d_coord_v2(smoothed_cl, delta);
        smoothed_cl = smoothed_cl_delta(:,2:4);
        
       
        centroids_LCCA = load(strcat(file_path,'/',patient_ID, mode_i,'_LCCA_centroid.txt'));
        %% Reference coordinate from LCCA origin
        % Reference coordinate of LCCA
        LCCA_reference_point = centroids_LCCA(1,:); 
        % Mark the location on the aortic innerline closest to the LCCA
        % marker.
        LCCA_innerline_marker = find_closest_point(smoothed_line,LCCA_reference_point);
        % Line for visual verification
        LCCA2inner_ref_line = reference_line(LCCA_reference_point,LCCA_innerline_marker);
        LCCA_centerline_marker = find_closest_point(smoothed_cl, LCCA_reference_point);

        %% Bird-beak algorithm
        if all_modes{k} == "post_ID"                       
            modelPlot = 4;               
            diameterPlot = [5 6];
            %% Load Graft Contours
            contours_graft = extract_contours(file_path,compose("seg_%s%s_endograft2",patient_ID, mode_i));
        
            %% Combine contours into one matrix
            combined_contours_graft = [];
            for n=1:length(contours_graft)
                %combine graft points into one matrix for easy plotting
                combined_contours_graft = [combined_contours_graft;contours_graft{n}];            
            end
            n_points_of_spacing = 1; % skip n points to conserve CPU
            graft_innerline = bruteforce_aortaline_contours(smoothed_line,contours_graft,n_points_of_spacing);
            
            g1 = graft_innerline(1,1:3); % Marks the inner point of the first contour
            a1 = find_closest_point(smoothed_line,g1);
            c1 = find_closest_point(smoothed_cl,g1);   
            %% Spline Inner Graft Line
        
            interpolated_graft_inner = cat(1,0,cumsum(sqrt(sum(diff(graft_innerline,[],1).^2,2))));
            graft_interp_inner = interp1(interpolated_graft_inner, graft_innerline, unique([interpolated_graft_inner(:)' linspace(0,interpolated_graft_inner(end),500)]),'PCHIP');
                
            % Endograft first contacts the aorta innerline
            
            origin_points = find_birdbeak_origin(smoothed_line, graft_interp_inner ,threshold);       
            a2 = origin_points(1,5:7);
            g2 = origin_points(1,2:4);        
            c2 = find_closest_point(smoothed_cl,a2);
            
            
            a1Fiducial = fiducialMarkerWrite(LCCA_innerline_marker, smoothed_line, a1);                      
            a2Fiducial = fiducialMarkerWrite(LCCA_innerline_marker, smoothed_line, a2);
                      
                     
            c1Fiducial = fiducialMarkerWrite(LCCA_centerline_marker, smoothed_cl, c1);                        
            c2Fiducial = fiducialMarkerWrite(LCCA_centerline_marker,smoothed_cl,c2); 
            
                      
            [bba, bbl] = birdbeak_length_angle2(g1,g2,a1,a2);
            
            bbh = norm(g1-a1); % bird-beak height
            bbh_all = [bbh_all;bbh];
            % Create line for visualization
            bbline = reference_line(g1,g2);
            output_LZ_preop(1,:) = [a1Fiducial, a2Fiducial, c1Fiducial, c2Fiducial];

        elseif all_modes{k} == "pre_ID"
            modelPlot = 1;
            diameterPlot = [2 3];           
            a1 = fiducialMarkerRead(output_LZ_preop(1,1),smoothed_line, LCCA_innerline_marker);
            a2 = fiducialMarkerRead(output_LZ_preop(1,2),smoothed_line,LCCA_innerline_marker);              
            c1 = fiducialMarkerRead(output_LZ_preop(1,3),smoothed_cl, LCCA_centerline_marker);
            c2 = fiducialMarkerRead(output_LZ_preop(1,4),smoothed_cl, LCCA_centerline_marker);
            
        end
        c_lcca_ref = find_closest_point(smoothed_cl,LCCA_reference_point);

        %% Angulation and curvature of the aorta
        
        window = 30; % in mm, for taking curvature and angulation
        
        [angle_a1, curve_a1] = angulation_curvature_at_point(smoothed_line,a1,window);
        [angle_c1, curve_c1] = angulation_curvature_at_point(smoothed_cl,c1,window);
        [angle_a2, curve_a2] = angulation_curvature_at_point(smoothed_line,a2, window);
        [angle_c2, curve_c2] = angulation_curvature_at_point(smoothed_cl,c2,window);
  %      [angle, curve] = angulation_curvature_at_point(smoothed_line,plp,window_min);
        boufi_angle = boufi_angulation(smoothed_line,a1,window);
        boufi_angles = [boufi_angles;boufi_angle];
        
        %[angles_all, curves_all] = vary_curvature_and_angulation(smoothed_line_delta,window,window,1)
        [~, curves_innerline, ~] = curvature_and_angulation_with_arclength(smoothed_line_delta, window, window, 1);
        [~, curves_centerline, ~] =  curvature_and_angulation_with_arclength(smoothed_cl_delta, window, window, 1);

        % Crop off the first few aortic contours 
        contours_crop = contours_aorta(6:end-1);
        centroids_crop = centroids_aorta(6:end-1,:);
        
        %% 
        [diameters, smoothed_markers,smoothed_arclength, interpolated] = ...
            profile_diameter_vs_arclength_cl(smoothed_cl_delta,contours_aorta);
        
        
        
        k_gauss = fspecial('gaussian',[1,300],50);
        filtered_diameters = conv(interpolated(:,2),k_gauss, "valid");
        filtered_arclength = round(conv(interpolated(:,1),k_gauss, "valid"),1);
        filtered_arc_diam = [filtered_arclength, filtered_diameters];
        
        landingZoneMarkerOnCenterline = find_index_of_row_in_array(c1,smoothed_cl);
        landingZoneMarkerOnCenterline = smoothed_cl_delta(landingZoneMarkerOnCenterline,1);
        
        % Display window where we are obtaining the diameter and curvature
        % for plotting
        ptchidx = (curves_centerline(:,1) >= landingZoneMarkerOnCenterline) & (curves_centerline(:,1) <= landingZoneMarkerOnCenterline+window);                                  % Area To Shade
        ptchidx_prox = landingZoneMarkerOnCenterline;
        ptchidx_dist = curves_centerline(ptchidx,1);
        ptchidx_dist = ptchidx_dist(end,1);

        %% mean curvature on cl over 30 mm window
        
        innerline_curve_at_lz_over_30mm = measure_over_window(curves_innerline,smoothed_line_delta,a1,30);
        centerline_curve_at_lz_over_30mm = measure_over_window(curves_centerline,smoothed_cl_delta,c1,30);
        diameter_over_30mm_window = measure_over_window(filtered_arc_diam,smoothed_cl_delta,c1,30);
        
        %% graft deployment
        %[percent_expanded,reported] = measure_proximal_endograft_size(contours_graft,graft_sizing,p); 
        %graft_expansion_index = [graft_expansion_index;percent_expanded, reported];
        %% Display Headers
        fprintf('\nGenerating plot data...')
        header1 = compose("Patient ID: \t %s%s",patient_ID,mode_i);
        header2 = compose("Birdbeak height: %f", round(bbh,2));
        header3 = compose("Birdbeak angle: %f", round(bba,2));

        %% Plot the Results 
 
        subplot(2,3,modelPlot)
        plot3(smoothed_line(:,1),smoothed_line(:,2),smoothed_line(:,3),"Color",[0.49 0.18 0.56],'MarkerSize', 4,'LineWidth',3),hold on
%         plot3(graft_inner(:,1),graft_inner(:,2),graft_inner(:,3),'.'), hold on
        %plot3(combined_contours_aorta(:,1),combined_contours_aorta(:,2),combined_contours_aorta(:,3),'.','Color',[0.49 0.18 0.56],'MarkerSize',3),hold on
        plot3(combined_contours_aorta(:,1),combined_contours_aorta(:,2),combined_contours_aorta(:,3),'.',"Color",[0 0.4470 0.7410],'MarkerSize',3),hold on
        plot3(LCCA_reference_point(1,1),LCCA_reference_point(1,2),LCCA_reference_point(1,3),'r*'), hold on
        plot3(smoothed_cl(:,1),smoothed_cl(:,2),smoothed_cl(:,3),"-","Color",[0.8500 0.3250 0.0980],'LineWidth',3)
        plot3(c1(:,1),c1(:,2),c1(:,3),"*","Color",[0 0.45 0])
        
        plot3(c_lcca_ref(:,1),c_lcca_ref(:,2),c_lcca_ref(:,3),"k*")
        if all_modes{k} == "post_ID"
            plot3(combined_contours_graft(:,1),combined_contours_graft(:,2),combined_contours_graft(:,3),'.','Color',[0.93 0.69 0.13], 'MarkerSize',3), hold on
            title(compose("Patient %d, Post ID Model",p),"FontSize",16)
            postop_geometric_output_data(patient_count,:) = [p innerline_curve_at_lz_over_30mm{1} centerline_curve_at_lz_over_30mm{1} diameter_over_30mm_window{1} curve_a1 curve_c1]
            birdbeak_and_graft(patient_count,:) = [p bbh bbl bba graft_sizing(p)]

        elseif all_modes{k} == "pre_ID"
            title(compose("Patient %d, Pre ID Model",p),"FontSize",16)
            preop_geometric_output_data(patient_count,:) = [p innerline_curve_at_lz_over_30mm{1} centerline_curve_at_lz_over_30mm{1} diameter_over_30mm_window{1} curve_a1 curve_c1]
        end
        %title({header1,header2,header3},'Interpreter', 'none')
        xlabel('x (mm)',"FontSize",16)
        ylabel('y (mm)',"FontSize",16)
        zlabel('z (mm)',"FontSize",16)
        axis 'equal', grid on
        hold off
        
        subplot(2,3,diameterPlot);
        yyaxis left
        plot(filtered_arclength,filtered_diameters, "-",'LineWidth',1), hold on
        
        % Mark the bird-beak zone on CL
        %a1_cl_marker = find_closest_point(smoothed_cl,g1);
        
        
        
        
        plot([landingZoneMarkerOnCenterline;landingZoneMarkerOnCenterline],[0;max(diameters)*1.25],'k--','LineWidth',1)
        plot(smoothed_arclength,diameters,"o")
        
        
        LCCA_innerline_index = find_index_of_row_in_array(LCCA_centerline_marker,smoothed_cl);
        LCCA_arclength = smoothed_cl_delta(LCCA_innerline_index,1);
        %LCCA_innerline_index = find_index_of_row_in_array(LCCA_innerline_marker, smoothed_line);
        %LCCA_arclength = smoothed_line_delta(LCCA_innerline_index,1);
        plot([LCCA_arclength;LCCA_arclength],[0;max(diameters)*1.25],'k--','LineWidth',1)
        xlabel("Arclength distance along centerline (mm)","FontSize",16)
        ylabel("Diameter (mm)","FontSize",16)
        h = text(landingZoneMarkerOnCenterline-4,max(diameters)/5,"Landing Zone");
        set(h,'Rotation',90);
        h = text(LCCA_arclength-4,max(diameters)/5,"LCCA");
        set(h,'Rotation',90);
        %title(compose("Patient %i",p),"FontSize",16)
        
        hold off
        yyaxis right
        plot(curves_centerline(:,1),curves_centerline(:,2),"-",'LineWidth',1)
        ylim([0 0.07]);
        xlim([0 300]);
        ylabel("Curvature (mm^{-1})","FontSize",16)
        title("Curvature and diameter along centerline arclength","FontSize",16)
        hold on
        t = curves_centerline(:,1);
        y = curves_centerline(:,2);
        patch([ptchidx_prox ptchidx_prox ptchidx_dist ptchidx_dist], [0 0.7 0.7 0], [0.6 0.4 0.9], 'FaceAlpha',0.20, 'EdgeColor','none')

        grid on
        %% Finished with patient, p.
        fprintf('finished\n'), toc % record time interval
        fprintf('\n\n\n')
        message1 = strcat('Birdbeak height: \t',num2str(bbh) ,' mm\n');
        message2 = strcat('Birdbeak length: \t',num2str(bbl) ,' mm\n');
        message3 = strcat('Birdbeak angle: \t', num2str(bba), ' degrees\n');
        fprintf([message1, message2, message3], 'Interpreter', 'none')                
        hold off
    end
    %% Save Output
    patient_count = patient_count+1;
end

%% For Printing Output to Screen.
fprintf('\nData format:\nCol1:\t\tBirdbeak height\nCol2-n: \tCurvatures at given boundaries\n')






