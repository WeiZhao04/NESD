function f_NESD_fMRI(fname,opts)
% Perform the Native Eigenspace Denoising (NESD) for rs-fMRI data
% fname:    (str):the full path of the input nii file
% opts:     (struct):basic parameter setting to run NESD

%%
mask = f_spm_load_nii(opts.mask_dir);
mask = logical(mask);
%% load and sliding-window the data
Header = spm_vol(fname);
temp = f_spm_load_nii(Header);
% masking
temp = f_masking(temp,mask);
[P,L] = size(temp);
winL = opts.win_len;
stepL = opts.step_len;
NoSeg = 1:stepL:(L-winL+1);
N = length(NoSeg);
Weight_rc = zeros(1,L);
INX = cell(1,N);
for k = 1:N
    if k < N
        INX{k} = NoSeg(k):NoSeg(k)+winL-1;
        Weight_rc(INX{k}) = Weight_rc(INX{k}) + 1;
    else
        INX{k}= NoSeg(k):L;
        Weight_rc(INX{k}) = Weight_rc(INX{k}) + 1;
    end
end
%% load the features
load([opts.fea_dir filesep 'RP.mat'],'rp24');
prin_signal = load([opts.fea_dir filesep 'PS.mat'],'PS');
seg_roi = f_spm_load_nii([opts.fea_dir filesep 'seg_roi.nii']);
seg_roi = f_masking(seg_roi,mask);
%%
deE.vec = cell(1,N);
deE.lat = cell(1,N);
deE.cof = cell(1,N);
mVal = mean(temp,2); 
for k_seg = 1:N
    disp(['svd decomposing, ' num2str(k_seg) '/' num2str(N) ', ...'])
    seg = temp(:,INX{k_seg});
    seg = detrend(seg',0)';
    seg(isnan(seg)) = 0;
    %% segmentary SVD
    [U,S,V] = svd(seg,'econ');  %  variance computed as diag(S).^2/size(seg,1);
    %% Temporal Features
    % 1) motion variation
    seg_RP = rp24(INX{k_seg},1:24);
    TFeature.RPcorr = corr(V,seg_RP);
    % 2) wm csf signal detection, aCompCorr
    seg_PS = prin_signal.PS(INX{k_seg},:);
    TFeature.PScorr = corr(V,seg_PS);
    % 3) ultra high/low powerspec
    [TFeature.spec,power_reg] = f_spectrum_metric(V,opts.TR);
    %% Spatial Features / not solid algorithm, uesd as qc currently
    % 1) gm,wm,csf,skull fraction
    [SFeature] = f_region_fraction(U,seg_roi,mask);
    %% detection and denoise
    [deV,otl_V] = f_regress_nv(V,TFeature,[seg_RP,seg_PS],power_reg);
    [deU] = f_regress_nu(U,seg_roi,mask);
    save(fullfile(opts.fea_dir,sprintf('outliers%03i.mat',k_seg)),'otl_V','TFeature','SFeature');
    deE.vec{k_seg} = deU;
    deE.lat{k_seg} = diag(S);%     E.vals{ks} = diag(S).^2/size(deAvg,1);
    deE.cof{k_seg} = deV;
    %% qality assesment metrics / turn off if not needed
    % 1) autocorr (pre/post)
    [QC.auto,QC.auto_post] = f_temporal_autocorr(V,deV);
    % 2) kurtosis,skewness,entropy (pre/post)
    [QC.VD,QC.VD_post] = f_volume_distribution(U,deU,mask);
    % 3) sepctrum (pre/post)
    QC.spec = TFeature.spec;
    QC.spec_post = f_spectrum_metric(deV,opts.TR);
    % 4) gm,wm,csf,skull fraction (pre/post)
    QC.frac = SFeature;
    [QC.frac_post] = f_region_fraction(deU,seg_roi,mask);
    save(fullfile(opts.fea_dir,sprintf('QC%03i.mat',k_seg)),'QC');
end
% free the memory
clear temp
clear U S V
%%
disp('Applying zero-block diag matrix for adjacency graph')
[rbook,TimePoint] = f_blkdiagMask(deE.vec);
%%
disp('graphic constructing')
[NH] = f_NHestablish_indiv(rbook,TimePoint);
%%
disp('GLS approximating')
[Y] = f_GSL_appr_indiv(NH,deE);
% opts.Len = NH.Len;
opts.kept_index = NH.KI;  % index of eigenvectors that kept by graphic adjacency matrices
%% 
disp('denoise segments reconstruction')
temp_full = zeros(P,L);
for k = 1:N
    RC = Y.vec{k}*diag(Y.lat{k})*Y.cof{k}'+mVal;
    temp_full(:,INX{k}) = temp_full(:,INX{k}) + RC;
end
temp_full = bsxfun(@rdivide,temp_full,Weight_rc);  % transfer to single to save memory
[dx,dy,dz] = size(mask);
Vol_new = single(zeros(L,dx,dy,dz));
Vol_new(:,mask>0) = temp_full';
Vol_new = permute(Vol_new,[2,3,4,1]);
outfile = [opts.out_dir filesep 'NESD_fMRI.nii'];
f_spm_save_nii(Vol_new,outfile,opts.mask_dir);
save(fullfile(opts.out_dir,'opts.mat'),'opts')
%%

end


%%%%%%%%%%%%%%%%%%%%% masking the data %%%%%%%%%%%%%%%%%%%%%%
function Mdata = f_masking(data,mask)

[Xm,Zm,Ym] = size(mask);

[Xd,Zd,Yd,Td] = size(data);

if isequal(Xm,Xd) && isequal(Zm,Zd) && isequal(Ym,Yd)
    data = reshape(data,[],Td);
    Mdata = data(mask>0,:);
else
    error('Dimension of mask is not matched with data %s',num2str([Xm,Zm,Ym]));
end

end

%%%%%%%%%%%%%%%%%%%%% THE NEW ONE (memory effiveint for computation cost) %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% generate a zero block diag matrix and mask with it %%%%%%%%%%%%%%%%%%%%%

function [rbook,TP] = f_blkdiagMask(Evecs)

n = size(Evecs{1},1);
N = size(Evecs,2);
rbook = cell(N,N);
TP = zeros(N,1);
for k = 1:N
    Evecs{k} = bsxfun(@minus,Evecs{k},sum(Evecs{k},1)/n);
    TP(k) = size(Evecs{k},2);
end

for kr = 1:N
    for kn = 1:N
        if kr <= kn
%             rbook{kr,kn} = [];
        else
            rbook{kr,kn} = Evecs{kr}'*Evecs{kn};
        end
    end
    disp(['calculating distance, ' num2str(kr) '/' num2str(N) ', ...'])
end

end

