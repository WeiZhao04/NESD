function f_plot_time_freq(fname,seg_roi,TR,winL,stepL)
% The creation of grayplot and relative FD and global signal lineplot.
% fname:    (str) - the path and filename of the fMRI file (*.img, *.hdr, *.nii, *.nii.gz)
% seg_roi:  (str) - the path of gm and wm ROI
% TR:       (int): the TR in seconds, to perform band-pass filter.
%           default 0.01-0.08Hz
% winL:     (int): the sliding window length, default 20 TRs(at least 20
%           volumes or 30 seconds)
% stepL:    (int): the step length, default 2 TRs(at least 2 TRs or 5 seconds)
%

[out_dir,name,ext,~] = spm_fileparts(fname);

mni_data = f_spm_load_nii(fname);
mni_data = mni_data.*any(mni_data,4).*any(seg_roi,4);
mni_temp = reshape(mni_data,[],size(mni_data,4));
clear mni_data

gm = mni_temp(seg_roi(:,:,:,1)>0,:);
wm = mni_temp(seg_roi(:,:,:,2)>0,:);
% csf = mni_temp(seg_roi(:,:,:,3)>0,:);
clear mni_temp
wb = cat(1,gm,wm);
wb(~any(wb,2),:) = [];
gl_signal = mean(wb,1,'omitnan');
gl_signal = y_IdealFilter(gl_signal',TR,[0.01 0]);
%% time-frequency
L = length(gl_signal);
if  L < 1024
    N = 1024;
else
    N = 2^nextpow2(L);
end
gl_signal = gl_signal/prctile(gl_signal,95);
[~,F_raw,T_raw,P_raw] = spectrogram(gl_signal,winL,winL-stepL,N,1/TR);
figure('visible','off'),
imagesc(T_raw,F_raw,P_raw),caxis([0 10])
set(gca,'yticklabel',[],'xticklabel',[],'looseInset',[0 0 0 0])
box on
set(gcf,'unit','centimeters','innerPosition',[5 5 4.5 2.4])
out_name = strcat('time_freq_',name,'.tiff');
print(gcf,[out_dir,filesep out_name], '-dtiff','-r600' );

end