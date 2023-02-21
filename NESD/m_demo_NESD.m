clc;clear;
close all;
%% readme before run
% (1) pls organize the data into one working folder first, and each subfolder
% better named as a particluar ordered patterns. Each subfolder only nedd
% to contain a signal native space file, named as fMRI.nii
% (2) Be sure that required code package, spm12 are added into PATH.
% e.g addpath(genpath('D:\code_packages\spm12'))
% (3) Be sure that all NESD codes are added into PATH and do not rename or
% delete 'm_demo_run.m'
% [codepath, ~, ~] = fileparts(which('m_demo_NESD.m'));
% e.g addpath(genpath(codepath))
%% 
wd = ['DATA']; % the working directionary
cd(wd)
sub_list = dir([wd filesep '*_*']); % the separator used to dected subfolders
%%
for k_sub = 1:length(sub_list)
    data_list = sub_list(k_sub).name;
    ori_fMRI = [wd filesep data_list filesep 'fMRI.nii'];
    [ori_dir,fun_name,ext,~] = spm_fileparts(ori_fMRI);
    %% parameters setting
    % setup the TR in seconds
    TR = 3;
    % setup the sliding-window parameters
    window_length = 50;
    step_length = 20;
    %% default subfolers creation
    % setup the preproc path
    pre_dir = [ori_dir filesep 'data_pre'];
    if ~exist(pre_dir,'dir')
        mkdir(pre_dir)
    else
        rmdir(pre_dir,'s')
        mkdir(pre_dir)
    end
    % setup the feature path
    fea_dir = [ori_dir filesep 'features'];
    if ~exist(fea_dir,'dir')
        mkdir(fea_dir)
    else
        rmdir(fea_dir,'s')
        mkdir(fea_dir)
    end
    % setup the output path
    out_folder = 'Result';
    out_dir = [ori_dir filesep out_folder];
    if ~exist(out_dir,'dir')
        mkdir(out_dir)
    else
        rmdir(out_dir,'s')
        mkdir(out_dir)
    end
    %% load raw data
    V = spm_vol(ori_fMRI);
    ori_data = f_spm_load_nii(ori_fMRI);
    [xd,yd,zd,td] = size(ori_data);
    %% ******************pre analysis**************************
    f_spm_preproc2(ori_fMRI,pre_dir);
    %% Region defination: mask generation for gray matter,white matter,csf,edge,sinus region
    seg_roi = zeros(xd,yd,zd,5);
    pro_roi = zeros(xd,yd,zd,5);
    prob = [0.2,0.2,0.2,0.2,0.9]; % a soft threhold for mask generation
    for k_roi = 1:5
        tmp_roi = f_spm_load_nii([pre_dir filesep 'c' num2str(k_roi) 'mean_' fun_name ext]);
        pro_roi(:,:,:,k_roi) = tmp_roi;
        seg_roi(:,:,:,k_roi) = logical(tmp_roi>prob(k_roi));
    end
    seg_roi_name = [fea_dir filesep 'seg_roi.nii'];
    f_spm_save_nii(seg_roi,seg_roi_name,ori_fMRI);
    %% generate brain mask in natvie space
    tempVol = logical(sum(seg_roi(:,:,:,1:3),4));
    mask = f_Automask(tempVol);
    mask_name = [ori_dir filesep 'raw_mask.nii'];
    f_spm_save_nii(mask,mask_name,ori_fMRI);
    clear tempVol
    %% aCompCorr, top 5 components
    ori_temp = reshape(ori_data,[],size(ori_data,4));
    ori_temp = ori_temp(mask>0,:);
    pro_temp = reshape(pro_roi,[],5);
    pro_temp = pro_temp(mask>0,:);
    wm = ori_temp(pro_temp(:,2)>0.9,:); % a relative high threshold for confidence
    csf = ori_temp(pro_temp(:,3)>0.9,:);
    [~,~,PC] = svd(cat(1,wm,csf),'econ');
    PS = PC(:,1:5);
    save(fullfile(fea_dir,'PS.mat'),'PS')
    %% computed GS and reigonal mean time-course in native space
    regional_signals = ori_temp'*pro_temp;
    regional_signals = bsxfun(@minus, regional_signals, mean(regional_signals));
    % global signal computed with the mean of gm, wm and csf
    regional_signals(:,6) = mean(regional_signals(:,1:3),2); 
    save(fullfile(fea_dir,'regional_signals.mat'),'regional_signals')
    save([fullfile(fea_dir,'gm_wm_csf_edge_sinus_GS_ordered'), '.txt'],'regional_signals', '-ascii', '-double','-tabs')
    clear ori_temp pro_temp wm_csf gm cr PC  % free memroy
    %% Frison 24 and Frame-wise displacement based on motion estimates
    rp = importdata(fullfile(ori_dir,'rp_fMRI.txt'));
    rp24 = [rp,[zeros(1,size(rp,2));rp(1:end-1,:)], rp.^2, [zeros(1,size(rp,2));rp(1:end-1,:)].^2];
    [rel_rms, abs_rms] = f_FD_Jenkinson(rp,V(1).fname);
    save(fullfile(fea_dir,'RP.mat'),'rp','rp24','abs_rms','rel_rms')
    %% (optional) computed original space SVD (!!! remove variable name in line 73)
%     [u,s,v] = svd(ori_temp,'econ'); 
%     temp = zeros(xd*yd*zd,td);
%     temp(mask>0,:) = zscore(u,0,2);
%     temp = reshape(temp,xd,yd,zd,td);
%     f_spm_save_nii(temp,fullfile(fea_dir,'Native_SVD_Umap.nii'),ori_fMRI);
%     save(fullfile(fea_dir,'Native_SVD_SV.mat'),'s','v')
%     clear temp u s v
    %% *************NESD section*****************
    opts.win_len = window_length;
    opts.step_len = step_length;
    opts.mask_dir = mask_name;
    opts.out_dir = out_dir;
    opts.fea_dir = fea_dir; 
    opts.TR = TR;
    f_NESD_fMRI(fullfile(ori_dir,['r' fun_name ext]),opts);
    nesd_name = fullfile(out_dir,strcat('NESD_',fun_name,ext));
    %% *****************post analysis****************************
    %% motion estimation for denoised data
    nesd_V = spm_vol(nesd_name);
    P = spm_realign(nesd_V,[]); % change the file if needed
    n = length(P);
    rp = zeros(n,6);
    for j=1:n
        qq = spm_imatrix(P(j).mat/P(1).mat);
        rp(j,:) = qq(1:6);
    end
    [rel_rms, abs_rms] = f_FD_Jenkinson(rp,nesd_V(1).fname);
    save(fullfile(out_dir,'NESD_RP.mat'),'rp','abs_rms','rel_rms')
    %% gray plot, better visulization effect for smoothed raw data
    fname = [ori_dir filesep 'fMRI.nii'];
    load([fea_dir filesep 'RP.mat'],'rel_rms')
    f_grayplot_rp_gs(fname,seg_roi,rel_rms);
    f_plot_time_freq(fname,seg_roi,TR,0.5*window_length,0.1*step_length);
    
    fname = [out_dir filesep 'NESD_fMRI.nii'];
    load([out_dir filesep 'NESD_RP.mat'],'rel_rms')
    f_grayplot_rp_gs(fname,seg_roi,rel_rms);
    f_plot_time_freq(fname,seg_roi,TR,0.5*window_length,0.1*step_length);
    %% (optional) Normalise with EPI, better performed with other standard protocals/pipelines
%     [codepath, ~, ~] = fileparts(which('m_demo_NESD.m'));
%     [SPMFilePath, ~, ~] = fileparts(which('spm.m'));
%     SPMJOB = load([codepath filesep 'templates' filesep 'Jobmats' filesep 'Normalize.mat']);
%     mean_fun = [pre_dir filesep 'mean_' fun_name ext];
%     FileList = [{ori_fMRI};{mean_fun};{seg_roi_name};{mask_name};{nesd_name}];
%     SPMJOB.matlabbatch{1,1}.spm.spatial.normalise.estwrite.subj(1,1).source={mean_fun};
%     SPMJOB.matlabbatch{1,1}.spm.spatial.normalise.estwrite.subj(1,1).resample=FileList;
%     SPMJOB.matlabbatch{1,1}.spm.spatial.normalise.estwrite.eoptions.template={[SPMFilePath,filesep,'toolbox',filesep,'OldNorm',filesep,'EPI.nii,1']};
%     SPMJOB.matlabbatch{1,1}.spm.spatial.normalise.estwrite.roptions.bb = [-90 -126 -72; 90 90 108];
%     SPMJOB.matlabbatch{1,1}.spm.spatial.normalise.estwrite.roptions.vox = [3 3 3];
%     spm_jobman('run',SPMJOB.matlabbatch);
    %% (optional) Smooth
%     SPMJOB = load([codepath filesep 'templates' filesep 'Jobmats' filesep 'Smooth.mat']);
%     w_raw = fullfile(ori_dir,strcat('w',fun_name,ext));
%     SPMJOB.matlabbatch{1,1}.spm.spatial.smooth.data = [{w_raw}];
%     SPMJOB.matlabbatch{1,1}.spm.spatial.smooth.fwhm = [6 6 6];
%     spm_jobman('run',SPMJOB.matlabbatch);
%     
%     fname = [out_dir filesep 'wNESD_fMRI.nii'];
%     [~,~,FWHM_est, ~] = y_Smoothest(fname,[ori_dir filesep 'wraw_mask.nii']);
%     FWHM = real(sqrt([9 9 9].^2-FWHM_est.^2));
%     SPMJOB = load([codepath,'templates',filesep,'Jobmats',filesep,'Smooth.mat']);
%     w_mni = fullfile(out_dir,strcat('wNESD_',fun_name,ext));
%     SPMJOB.matlabbatch{1,1}.spm.spatial.smooth.data = [{w_mni}];
%     SPMJOB.matlabbatch{1,1}.spm.spatial.smooth.fwhm = FWHM;
%     spm_jobman('run',SPMJOB.matlabbatch);
    %% (optional) static and dynamic FC, only availalbe in MNI space (2mm or 3mm)
%     altas_name = [codepath filesep 'templates' filesep 'altas' filesep Schaefer2018_100Parcels_7Networks_order_FSLMNI152_3mm.nii'];
%     mni_seg_roi = f_spm_load_nii([fea_dir filesep 'wseg_roi.nii']); 
%     fname = [ori_dir filesep 'swfMRI.nii'];
%     load([fea_dir filesep 'RP.mat'],'rel_rms')
%     f_grayplot_rp_gs(fname,mni_seg_roi,rel_rms);
%     f_plot_FC(fname,altas_name,TR,0.5*window_length,0.1*step_length);
%     f_plot_time_freq(fname,mni_seg_roi,TR,0.5*window_length,0.1*step_length);
% 
%     fname = [out_dir filesep 'swNESD_fMRI.nii'];
%     load([out_dir filesep 'NESD_RP.mat'],'rel_rms')
%     f_grayplot_rp_gs(fname,mni_seg_roi,rel_rms);
%     f_plot_FC(fname,altas_name,TR,0.5*window_length,0.1*step_length);
%     f_plot_time_freq(fname,mni_seg_roi,TR,0.5*window_length,0.1*step_length);    
end


