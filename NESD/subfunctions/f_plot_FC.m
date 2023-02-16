function [ROICorrelation,FCD] = f_plot_FC(fname,altas_name,TR,winL,stepL)
% Compute the staic FC(sFC) and FC dynamic(FCD)
% f_plot_FC(fname,altas_name,TR,winL,stepL)
% fname:    (str):the full path of the input nii file
% altas:    (str):the full path of the altas file of ROIs, it must be the
%           Schaefer's 7 networks(100 ROI version), or else change the plot
% TR:       (int): the TR in seconds, to perform band-pass filter.
%           default 0.01-0.08Hz
% winL:     (int): the sliding window length, default 20 TRs(at least 20
%           volumes or 30 seconds)
% stepL:    (int): the step length, default 2 TRs(at least 2 TRs or 5 seconds)
%
% optional(not code yet)


%%  compute the static FC
[out_dir,name,~,~] = spm_fileparts(fname);
ROI_altas = f_spm_load_nii(altas_name);

mni_data = f_spm_load_nii(fname);
TP = size(mni_data,4);
mni_temp = reshape(mni_data,[],TP);
clear mni_data
mni_temp = mni_temp(ROI_altas>0,:);
mni_temp(isnan(mni_temp)) = 0;
% detrend
mni_temp = detrend(mni_temp',1)';

% filter 
mni_temp = y_IdealFilter(mni_temp',TR,[0.01 0.08]);

% FC
ROI_altas = reshape(ROI_altas,1,[]);
ROI_altas = ROI_altas(ROI_altas>0);

Label = unique(ROI_altas);
Label(isnan(Label)) = []; 
Label(Label==0) = []; 
RN = length(Label);
ROI_mCourses = zeros(TP,RN);
for k_i = 1:length(Label)
    ROI_mCourses(:,k_i) = mean(mni_temp(:,ROI_altas==Label(k_i)),2);
end
ROISignals = double(ROI_mCourses);
clear mni_temp
%% print sFC map
figure('visible','off'),
imagesc(corr(ROISignals)),caxis([-1 1])
hold on
h = plot([9.5,9.5],[0.5,100.5],[15.5,15.5],[0.5 100.5],[23.5,23.5],[0.5 100.5], ...
    [30.5,30.5],[0.5 100.5],[33.5,33.5],[0.5 100.5],[37.5 37.5],[0.5 100.5], ...
    [50.5,50.5],[0.5 100.5]);
set(h,'Color','k','linewidth',0.7)
hold on
h = plot([0.5,100.5],[9.5,9.5],[0.5 100.5],[15.5,15.5],[0.5 100.5],[23.5,23.5], ...
    [0.5 100.5],[30.5,30.5],[0.5 100.5],[33.5,33.5],[0.5 100.5],[37.5 37.5], ...
    [0.5 100.5],[50.5,50.5]);
set(h,'Color','k','linewidth',0.7)
% xticks([5 12 18 26 31.5 36.5 44 75])
% xticklabels({'Vis','SM','DA','SA','LT','CN','DMN','Right Sphere'})
set(gca,'ytick',[],'xtick',[],'looseInset',[0 0 0 0])
set(gcf,'unit','centimeters','innerPosition',[5 5 4.6 4.6])
axis off
out_name = strcat('FC_plot_',name,'.tiff');
print(gcf,[out_dir filesep out_name], '-dtiff','-r600' );
%% save ROISignals
save([fullfile(out_dir,['ROISignals_', name]), '.txt'], 'ROISignals', '-ASCII', '-DOUBLE','-TABS')

ROICorrelation = corrcoef(ROISignals);
save([fullfile(out_dir,['ROICorrelation_', name]), '.txt'], 'ROICorrelation', '-ASCII', '-DOUBLE','-TABS')

ROICorrelation_FisherZ = 0.5 * log((1 + ROICorrelation)./(1- ROICorrelation));
save([fullfile(out_dir,['ROICorrelation_FisherZ_', name]), '.txt'], 'ROICorrelation_FisherZ', '-ASCII', '-DOUBLE','-TABS')

%% compute the FC dynamic
mask = tril(ones(RN,RN),-1);
NoSeg = 1:stepL:(TP-winL+1);
N = length(NoSeg);
for k_seg = 1:N
    if k_seg < N
        FC_t(k_seg,:,:) = corr(ROISignals(NoSeg(k_seg):NoSeg(k_seg)+winL-1,:));
    else
        FC_t(k_seg,:,:) = corr(ROISignals(NoSeg(k_seg):TP,:));
    end
end
FCD = corr(FC_t(:,mask>0)');
%%
figure('visible','off'),
imagesc(FCD,[-1 1]),colormap('jet')
set(gca,'ytick',[],'xtick',[],'looseInset',[0 0 0 0])
set(gcf,'unit','centimeters','innerPosition',[5 5 4.6 4.6])
axis off
out_name = strcat('FCD_',name,'.tiff');
print(gcf,[out_dir filesep out_name], '-dtiff','-r600' );

end


