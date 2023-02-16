function [deU] = f_regress_nu(U,seg_roi,mask)
%  To perform the denoised of spatial eigenvectors by excluding small
% clusters and soft smoothing, so they are shaped into a better
% distribusion in revealing meaningful signals.
% U:   (matrice) the iuput spatial eigenvectors, should be reformed into 
%       voxel x samples
% seg_roi:(matrice) the 4D binary value matrice of gray matter, white
%       matters and csf; basicly scatter located in csf are rejected
% mask: (matrice) a 3D binary value matrice of whole barin mask
%  Outputs
% deU:  (matrice) the output denoised spatial eigenvectors

%%
% [dL ds] = size(U);
ConnectivityCriterion = 18;
ClusterSize = 50;
[P,N] = size(U);
for k_v = 1:N
    eigenvector = U(:,k_v);
    % transfer voxels in U vectors into volume
    Vol = zeros(size(mask));
    Vol(mask) = zscore(eigenvector);
    %         % smooth
    Vol = smooth3(Vol,'gaussian',[3 3 3]);
    % threshold and search clusters
    Vol_th = Vol .* ((Vol <= -1.65) + (Vol >= 1.65));
    [CMask, CNum] = bwlabeln(Vol_th,ConnectivityCriterion);
    %%
    for k_label = 1:CNum
        CLabel = CMask==k_label;
        numVoxels =  length(find(CLabel));
        peak_map = zeros(size(mask));
        % remove small clusters
        if numVoxels < ClusterSize
            eigenvector(CLabel(mask)) = 0;
        else % polar point detection
            Cintensity = Vol(CLabel);
            pos = find(CMask==k_label);
            peakpos = find(abs(Cintensity) == max(abs(Cintensity)));
            peakintensity = Cintensity(peakpos);
            [i,j,k] = ind2sub(size(Vol),pos(peakpos));
            peakind = [i,j,k];
            peak_map(i,j,k) = 1;
            peak_fraction = seg_roi'*peak_map(mask);
            cluster_frac = (seg_roi'*CLabel(mask)/length(pos));
            if peak_fraction(3)*cluster_frac(3)
                eigenvector(CLabel(mask)) = 0;
            end
        end
    end
%     [~,~,r1,~] = regress(eigenvector,reg_eigen);
    Vol = zeros(size(mask));
    Vol(mask) = eigenvector;
    % soft smooth
    eig_th = 0.15;
    neg_th = norminv(eig_th,mean(eigenvector),std(eigenvector));
    pos_th = norminv(1-eig_th,mean(eigenvector),std(eigenvector));
    Vol_soft = Vol.*((Vol <= neg_th) + (Vol >= pos_th));
    soft_norm = norm(Vol_soft(:));
    Vol_soft = smooth3(Vol_soft,'gaussian',[3 3 3]);
    Vol_soft = Vol_soft/norm(Vol_soft(:))*soft_norm;
    Vol(Vol_soft~=0) = Vol_soft(Vol_soft~=0);
    % norm the eigenvectors
    deU(:,k_v) = col_norm(Vol(mask));
end


%     eig_mask = ((Vol <= eig_th) + (Vol >= -eig_th));
%     Vol_soft = Vol .* eig_mask;
%     Vol_soft = smooth3(Vol_soft,'gaussian',[3 3 3]);
%     Vol(eig_mask>0) = Vol_soft(eig_mask>0);

%%
function X_norm = col_norm(X)
        [r,c] = size(X);
        X_norm = zeros(r,c);
        for col_num = 1:c
            temp = X(:,col_num);
            temp = temp - mean(temp);
            temp = temp/norm(temp);
            X_norm(:,col_num) = temp;
        end

end

end

