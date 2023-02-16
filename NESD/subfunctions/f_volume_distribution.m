function [VD,VD_post] = f_volume_distribution(U,deU,mask)
%  To compute the spatial eigenvectors' distribusion features including  
% Kurtosis, skewness, entropy
% U:   (matrice) the iuput spatial eigenvectors, should be reformed into 
%       voxel x samples
% deU: (matrice) the denoised spatial eigenvectors, should be reformed into 
%       voxel x samples
% mask: (matrice) a 3D binary value matrice of whole barin mask
%%
[~,NoV] = size(U);

for k_v = 1:NoV
    % transfer voxels in U vectors into volume
    Vol = zeros(size(mask));
    Vol(mask) = zscore(U(:,k_v));
    % the pre Kurtosis, skewness, entropy, variance
    VD.kur(k_v) = kurtosis(Vol(:));
    VD.ske(k_v) = skewness(Vol(:));
    VD.ent(k_v) = entropy(Vol(:));
    % the post
    Vol = zeros(size(mask));
    Vol(mask) = zscore(deU(:,k_v));
    %
    VD_post.kur(k_v) = kurtosis(Vol(:));
    VD_post.ske(k_v) = skewness(Vol(:));
    VD_post.ent(k_v) = entropy(Vol(:));
    
%     VM.KR(k_v) = abs(VM.postK-VM.preK)/VM.preK;
%     VM.SR(k_v) = abs(VM.postS-VM.preS)/VM.preS;
%     VM.ER(k_v) = abs(VM.postE-VM.preE)/VM.preE;
end
