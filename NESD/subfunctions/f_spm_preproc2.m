function f_spm_preproc2(fun_file,pre_dir)
% The second edition of segmentation on fMRI(w/o T1 assistance) based on spm12.
% fun_file: (str) - the path and filename of the fMRI file (*.img, *.hdr, *.nii, *.nii.gz)
% pre_dir:  (str) - the path to save output files for preproc

%%
[codepath, ~, ~] = fileparts(which('m_demo_run.m'));
[SPMPath, ~, ~] = fileparts(which('spm.m'));
[fMRI_dir,fMRI_name,ext,~] = spm_fileparts(fun_file);
%% 0) mean reference
[raw_fMRI] = f_spm_load_nii(fun_file);
mean_fMRI = mean(raw_fMRI,4);
mean_funfile = fullfile(fMRI_dir,strcat('mean_',fMRI_name,ext));
f_spm_save_nii(mean_fMRI,mean_funfile,fun_file);
copyfile(mean_funfile,pre_dir)
%% 1) Realignment
SPMJOB = load([codepath,'templates',filesep,'Jobmats',filesep,'Realign.mat']);
SPMJOB.matlabbatch{1,1}.spm.spatial.realign.estwrite.data{1,1}={fun_file};
spm_jobman('run',SPMJOB.matlabbatch);
delete(fullfile(fMRI_dir,[fMRI_name '.mat']))
%% 2) Dartel Segment

SPMJOB = load([codepath,'templates',filesep,'Jobmats',filesep,'NewSegment.mat']);
SPMJOB.matlabbatch{1,1}.spm.tools.preproc8.channel.vols={[pre_dir filesep 'mean_' fMRI_name ext]};
for tclass =1:6 
SPMJOB.matlabbatch{1,1}.spm.tools.preproc8.tissue(1,tclass).tpm{1,1}=[SPMPath,filesep,'tpm',filesep,'TPM.nii',',',num2str(tclass)]; %YAN Chao-Gan, 161006.
SPMJOB.matlabbatch{1,1}.spm.tools.preproc8.tissue(1,tclass).warped = [0 0]; % Do not need warped results. Warp by DARTEL
end
SPMJOB.matlabbatch{1,1}.spm.tools.preproc8.warp.affreg = 'mni';
preproc = SPMJOB.matlabbatch{1,1}.spm.tools.preproc8;
preproc.warp.mrf = 1;
preproc.warp.cleanup = 1;
preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
preproc.warp.fwhm = 0;
SPMJOB=[];
SPMJOB.matlabbatch{1,1}.spm.spatial.preproc = preproc;
spm_jobman('run',SPMJOB.matlabbatch);
%%
end