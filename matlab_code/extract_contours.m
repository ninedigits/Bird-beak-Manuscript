function [f_contours] = extract_contours(file_path,contours_file)
%{

This function extracts the contour coordinates from SimVascular contour
files

INPUTS:     SimVascular contours file
OUTPUTS:    Cell array with individual contour files   


EXAMPLE INPUT: 

file = seg_T1pre_ID_aorta2

file contents:

/group/aorta2/47
47
num_smooth_modes 12 xhat {0.000000 0.856293 -0.516490} nrm {0.657364 0.389213 0.645280} num_subsample_pts 330 creation_method hand_drawn pathId 10000 pos {-94.853700 131.454163 -156.968131} creation_type levelset posId 47 num_hand_picked_pts 20
-105.105637 120.147118 -139.704132
-104.866196 119.863892 -139.777237
-104.620766 119.594147 -139.864563
-104.369942 119.338417 -139.965820
...
-105.563965 120.750755 -139.601318
-105.338448 120.443054 -139.645477

/group/aorta2/58
58
num_smooth_modes 12 xhat {0.000000 0.866310 -0.499506} nrm {0.652867 0.378362 0.656206} num_subsample_pts 301 creation_method hand_drawn pathId 10000 pos {-92.421758 132.879270 -154.553713} creation_type levelset posId 58 num_hand_picked_pts 18
-104.238686 125.123306 -138.324905
-104.108978 124.783447 -138.257980
...

EXAMPLE OUTPUT

f_contours = 1x55 cell array

f_contours{1} = 

-105.105637 120.147118 -139.704132
-104.866196 119.863892 -139.777237
-104.620766 119.594147 -139.864563
-104.369942 119.338417 -139.965820
...
-105.563965 120.750755 -139.601318
-105.338448 120.443054 -139.645477

%}


% close all
% clear all
%  
% contours_file = 'seg_T1pre_ID_aorta2'; 
% 
% 
% 
% file_path = '/Users/maxfrohlich/Desktop/stanford_code/MATLAB_development_2017/T1'; 



message = strcat('extracting contours: \t', contours_file, '......');
fprintf(message)
fid = fopen(file_path + contours_file);

all = {};
contours = [];
while feof(fid) == 0 % loop continues until reaching end of file
    
    tline = fgets(fid);
    try
        
        num = str2num(tline);   % tries to convert lines into numerical array
        if isequal(size(num),[1 3]) % only add the 1 x 3 matrices
            contours =[contours;num];
        else
            contours = [contours;[0 0 0]]; % adds rows of zeros when reaching next contour to separate contour points
        end
    catch

    end

end

%contours

% length(contours);
% 
indices = find(ismember(contours,[0, 0, 0], 'rows')); %finds index where rows of [0, 0, 0] appear

f_contours = {};
contour_indices = [];
for i=1:length(indices)-1

    try
        if indices(i)~=(indices(i+1)-1) %removes any duplicate [0 0 0] so that indices start and end at the beginning and end of contour
            contour_indices = [contour_indices; (indices(i))];  
        else
            continue
        end
    catch
    end
end
contour_indices;


for j=1:length(contour_indices)
    count = 0;
    start = contour_indices(j)+1;
    line = start;
    while contours(line,:) ~= [0 0 0]
        line = line + 1;
    end
    stop = line-1;
    f_contours = [f_contours,contours((start:stop),:)];
end

fclose(fid);
fprintf('extracted \n')
end




