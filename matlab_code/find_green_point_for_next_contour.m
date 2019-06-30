function[G_next, D_next] = find_green_point_for_next_contour(index,direction,A,kontur,G1)
%Function finding the Greenwich point for the next contour along the set direction
% Input:                index       (the current aortic cross section plane)
%                       direction   (With or against the blood flow, i.e. direction = 1: distal and direction = -1: proximal);
%                       A           (A Nx3 matrix, where N is the number of segmented aortic sections and the coordinates
%                                    stands for the position of the centroids in those sections);
%                       kontur      (the structure of all countor points for all sections);
%                       G1          (the current Greenwhich point).
% Output:               G_next (the Greenwich point on the next aortic section)
%                       D_next (
% Supporting scritps:   define_plane_normal_from_contour
%                       find_intersecting_vector
% Used in:              slim_analysis
%
% Author: Torbj?rn Lundh
% Dept. Mathematical Sciences, Chalmers university of Technology
% email: torbj?rn.lundh@chalmers.se
% 2016; Last revision: 02-Aug-2017
%------------- <<>> --------------
index_next = index+direction; %next index, a step distally or proximally
kkont_next=kontur(index_next).field; % the next contour in aorta
cA = A(index,:); %centerpoint of current section
cA_next=A(index_next,:); %the next counter center
N_G1 = cross((cA-G1),(cA - cA_next));  %compute vector perpendicular to (cA-G1) and (cA - cA_next)
N_G1 = N_G1/norm(N_G1); %normalize
NA_next = define_plane_normal_from_contour(kkont_next,cA_next); %compute normalized normal from next contour

D_next = find_intersecting_vector(N_G1, cA_next, NA_next, cA_next);
D_next=sign(D_next*(-cA+G1)')*D_next; %to make sure that the direction is close to the previous Green point

dot_product_with_D=zeros(1,length(kkont_next));
for j=1:length(kkont_next)
    Rad_vector = kkont_next(j,:) - cA_next; % the direction of all counter points seen from the center point
    Rad_vector = Rad_vector/norm(Rad_vector); %normalize
    dot_product_with_D(j) = D_next*Rad_vector';
end
[~,j_max] = max(dot_product_with_D); %pick the closest counter point
G_next = kkont_next(j_max,:); %Greenwich point taken
end