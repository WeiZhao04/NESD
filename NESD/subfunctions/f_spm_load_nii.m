function [Y,V] = f_spm_load_nii(P,mask)
% FORMAT [Y,XYZ] = spm_read_vols(P,mask)
% P - a char or cell array of filenames
% mask - implicit zero mask?
% V - a structure array containing image volume information
%__________________________________________________________________________


if nargin<2, mask = 0; end
if nargin<1, error('insufficient arguments'), end

V = spm_vol(P);

Y = spm_read_vols(V,mask);

end