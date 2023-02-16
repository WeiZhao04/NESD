function [SFeature] = f_region_fraction(U,seg_roi,mask)
%  To compute the spatial eigenvectors' fration that those valid voxel 
% located in which part of the brain.
% U:   (matrice) the iuput spatial eigenvectors, should be reformed into 
%       voxel x samples
% seg_roi:(matrice) the 4D binary value matrice of gray matter, white
%       matters and csf
% mask: (matrice) a 3D binary value matrice of whole barin mask
%%
[~,nVol] = size(U);
nClass = size(seg_roi,2);
volume_frac = zeros(nVol,nClass);
for k_v = 1:nVol
    % transfer voxels in U vectors into volume
    Vol = zeros(size(mask));
    zU = zscore(U(:,k_v));
    Vol(mask) = zU;
    % threshold and search clusters
    Vol = Vol .* ((Vol <= -1) + (Vol >= 1));
    % the thresholded whole brain voxel region location, because of the
    % overlap the total fraction is not solid 1
    total_voxel = sum(Vol~=0,'all');
    volume_frac(k_v,:) = seg_roi'*logical(Vol(mask))/total_voxel;
end
%%
SFeature.VF = volume_frac;



