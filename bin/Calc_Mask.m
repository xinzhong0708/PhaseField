function [mask] = Calc_Mask(phi,thickness)
%Threshold
thres = 1e-10;
mask  = 0*phi;
for ip = 1:size(phi,3)
    %Find the region
    % 1. Threshold to get core region
    mask_core     =  phi(:,:,ip) > thres;

    % 2. Create a structuring element (disk shape is good for isotropy)
    se            = strel('disk', thickness);

    % 3. Dilate the mask to include a rim around
    mask(:,:,ip)  = imdilate(mask_core, se);
end

end