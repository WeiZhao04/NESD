function [auto,auto_post] = f_temporal_autocorr(V,deV)
%  To compute the temporal eigenvectors' autocorrelation features including  
% lagOne, lagTwo, lagThree
% V:   (matrice) the iuput temporal eigenvectors, should be reformed into 
%       timepoints x samples
% deV: (matrice) the denoised temporal eigenvectors, should be reformed into 
%       timepoints x samples
%%
[L,N] = size(V); % length and number of timecourse

for k = 1:N
[cof,lags] = xcorr(V(:,k));

IndOne = find(lags==1);
auto.lagOne(k) = abs(cof(IndOne));
auto.lagTwo(k) = abs(cof(IndOne+1));
auto.lagThree(k) = abs(cof(IndOne+2));
end

for k = 1:N
[cof,lags] = xcorr(deV(:,k));

IndOne = find(lags==1);
auto_post.lagOne(k) = abs(cof(IndOne));
auto_post.lagTwo(k) = abs(cof(IndOne+1));
auto_post.lagThree(k) = abs(cof(IndOne+2));
end



end