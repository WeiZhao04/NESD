function Vn = f_spm_save_nii(Y,name,tmp)
% Write an image volume to disk, setting scales and offsets as appropriate
% FORMAT V = spm_write_vol(V,Y)
% V (input)  - a structure containing image volume information (see spm_vol)
% Y          - a one, two or three dimensional matrix containing the image voxels
% V (output) - data structure after modification for writing.
%
% Note that if there is no 'pinfo' field, then SPM will figure out the
% max and min values from the data and use these to automatically determine
% scalefactors.  If 'pinfo' exists, then the scalefactor in this is used.
%__________________________________________________________________________

if nargin<3, error('insufficient arguments'), end

ext = spm_file(name,'ext');
prec   = {'uint8','uint8','int16','int32','float32','float64','int8','uint16','uint32'};
conv   = {'logical','uint8','int16','int32','single','double','int8','uint16','uint32'};
sel = find(strcmpi(conv,deblank(class(Y))));  

switch ext
    case 'nii'
        Vn = spm_vol(tmp);
        Vn = Vn(1);
        if ndims(Y) < 4
            Vn.fname = name;
            Vn.dim = size(Y);
            Vn.dat = Y;
            Vn.dt(1) = spm_type(prec{sel});
            V = spm_write_vol(Vn,Y);
        elseif ndims(Y) == 4
            TP = size(Y,4);
            for k = 1:TP
                Vn.fname = name;
                Y_vol = squeeze(Y(:,:,:,k));
                Vn.dim = size(Y_vol);
                Vn.dat = Y_vol;
                Vn.dt(1)  = spm_type(prec{sel});
                Vn.n = [k,1];
                V(k) = spm_write_vol(Vn,Y_vol);
            end
        end
    case 'gz'
        name = name(1:end-3);
        Vn = spm_vol(tmp);
        Vn = Vn(1);
        if ndims(Y) < 4
            Vn.fname = name;
            Vn.dim = size(Y);
            Vn.dat = Y;
            Vn.dt(1) = spm_type(prec{sel});
            V = spm_write_vol(Vn,Y);
        elseif ndims(Y) == 4
            TP = size(Y,4);
            for k = 1:TP
                Vn.fname = name;
                Y_vol = squeeze(Y(:,:,:,k));
                Vn.dim = size(Y_vol);
                Vn.dat = Y_vol;
                Vn.dt(1) = spm_type(prec{sel});
                Vn.n = [k,1];
                V(k) = spm_write_vol(Vn,Y_vol);
            end
            
        end
        gzip(V(1).fname);
        delete(V(1).fname);
    otherwise
        error(['".' ext '" is not a recognised extension.']);
end




end