function [deV,outliers] = f_regress_nv(V,TF,regs_temporal,regs_spectrum)
%  To perform the denoised of temporal eigenvectors by regressing out
%  Friston 24 head motion parameters, WM&CSF signals(aCompCorr) and
%  ultra-low sepctrum trend.
% V:   (matrice) the iuput temporal eigenvectors, should be reformed into 
%       timepoints x samples
% TF:  (struct) the tempraol features worked as flags to identify noisy
% eigenvectors for denoising
% regs_temporal: (matrice) the temporal noise source 
% regs_spectrum: (matrice) the spectrum noise source 
%
%  Outputs
% deV:  (matrice) the output denoised temporal eigenvectors, should be reformed into 
%       timepoints x samples
% outliers:(struct) outliers that detected by thresholding or data-driven outliers method  

%%
% [ds,dl] = size(V);
%%
outliers.RP = abs(TF.RPcorr)>0.05;
outliers.PS = abs(TF.PScorr)>0.05;
% varied a lot in different segmentations or subjects, swith to threshold,
% can be improved into more robust way
% outliers.RP = reshape(isoutlier(TF.RPcorr(:),'gesd'),size(TF.RPcorr)); 
% outliers.PS = reshape(isoutlier(TF.PScorr(:),'gesd'),size(TF.PScorr));
% spec is much more reliable and can be fixed in latter filter 
outliers.SpecUL = isoutlier(TF.spec.UL,'gesd')';
% when add new regressors, be sure to sync regs_spectrum
%%
% construct matrics of outliers and regressors
O_temporal = cat(2,outliers.RP,outliers.PS);
% regs_temporal = detrend(regs_temporal,0);
O_spectrum = cat(2,outliers.SpecUL); % outliers.SpecfALFF,,outliers.SpecUH,outliers.SpecCR,outliers.SpecCO2,outliers.SpecVA
% compute the flag to identify eigenvectors needed for denoise
flag_TC = logical(sum(O_temporal,2));
flag_Spec = logical(sum(O_spectrum,2));
%%
[L,N] = size(V);
for k = 1:N
    y = V(:,k);
    regs = [];
    if flag_TC(k) 
        idx_t = O_temporal(k,:);
        regs = [regs regs_temporal(:,idx_t)];
    end
    if flag_Spec(k)
        idx_s = O_spectrum(k,:);
        fft_temp = regs_spectrum(:,k,idx_s);
        fft_temp = sum(fft_temp,2);
        regs_temp = real(ifft(fft_temp));
        regs_sbase = regs_temp(1:N,:); % not needed indeed
        
        regs = [regs regs_sbase];
    end
    if ~isempty(regs)
        [~,~,r1,~] = regress(y,regs);
        deV(:,k) = (r1);
    else 
        deV(:,k) = y;
    end
    
    
end


end