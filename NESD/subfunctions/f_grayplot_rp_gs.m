function f_grayplot_rp_gs(fname,seg_roi,rel_rms)
% The creation of grayplot and relative FD and global signal lineplot.
% fname: (str) - the path and filename of the fMRI file (*.img, *.hdr, *.nii, *.nii.gz)
% seg_roi:  (str) - the path of gm and wm ROI
% rel_rms:  (vec) -  the relative FD course

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
mean_wb = mean(wb,2,'omitnan');

wb = detrend(wb',0)';
wb = detrend(wb',1)';
grayplot = bsxfun(@rdivide,wb,mean_wb)*100;
grayplot = bsxfun(@times,grayplot,any(grayplot,2));
gl_signal = mean(wb,1,'omitnan');

dd = rms(diff(wb,1,2),1);
dvars = [mean(dd) dd];
clear wb


%%  lineplot of GS, DVARS, rFD
figure('visible','off'),
imagesc(grayplot,[-2 2]),colormap(gray)
set(gca,'ytick',[],'xtick',[],'looseInset',[0 0 0 0])
set(gcf,'unit','centimeters','innerPosition',[10 10 4.6 1.2])
axis off
out_name = strcat('grayplot_',name,'.tiff');
print(gcf,[out_dir filesep out_name], '-dtiff','-r600' );
%% grayplot
figure('visible','off'),
gl_signal = detrend(gl_signal);
gl_signal = gl_signal/max(abs(gl_signal))/2+1;
plot(gl_signal,'b','linewidth',1);
hold on
plot(rel_rms,'r','linewidth',1);
dvars = detrend(dvars);
dvars = dvars/max(abs(dvars))/2+1;
plot(dvars,'k:','linewidth',1);
xlim([0 length(gl_signal)])
ylim([0 2])
set(gca,'ytick',[],'xtick',[],'looseInset',[0 0 0 0])
set(gcf,'unit','centimeters','innerPosition',[10 10 4.6 1.5])
axis off
out_name = strcat('GS_rFD_DVARS_plot_',name,'.tiff');
print(gcf,[out_dir filesep out_name], '-dtiff','-r600' );

end