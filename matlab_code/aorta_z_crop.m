function aorta_zc = aorta_z_crop(aorta_s)
%this function crops the aorta, the centroids, and the contours 
start = aorta_s(1,3);


for i=1:length(aorta_s)
    if aorta_s(i,3)<start
        stop = i-1;
        break
    else
        stop = i-1;
    end
end

aorta_zc = aorta_s(1:stop,:);

