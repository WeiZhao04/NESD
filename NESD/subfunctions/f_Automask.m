function [mask Header_Out] = f_Automask_t1(Data)

% function [mask Header_Out] = w_Automask(Data, OutputFile, Verb)
% Input:
% Data - Could be the following format:
%                  1. A single file (.img/hdr, .nii, or .nii.gz), give the path and filename.
%                  2. A directory, under which could be a single 4D file, or a set of 3D images
%                  3. A Cell (nFile * 1 cells) of filenames of 3D image file, or a single file of 4D NIfTI file.
% OutputFile - The output file name of the mask.
% Verb - 1: output the details; 0: no output of details.
% Output: the Automasks.

%%
% Revised by Xin-di Wang, 140106. 
% State Key Laboratory of Cognitive Neuroscience and Learnling,
% Beijing Normal University, Beijing, PR China
% sandywang.rest@gmail.com
%                                    
% Modified according to AFNI's 3dAutomask.
%%
mask = Data;



[nx , ny , nz]=size(mask);

%Select the Largest Cluster
mask = w_mask_clust(mask);




%Select Largest Cluster again

mask = w_mask_clust(mask);

%Erode
npeel=1;
peelthr=17;
mask=w_mask_erodemany(mask, npeel, peelthr);


%Fill small holes
[mask , ii]=w_mask_fillin_once(mask , 1);
jj=ii;
if ii>0
    [mask , ii]=w_mask_fillin_once(mask ,1);
    jj=jj+ii;
    if ii>0
        [mask , ii]=w_mask_fillin_once(mask , 1);
        jj=jj+ii;
    end
end

mask=w_mask_erodemany( mask, npeel, peelthr);

mask = w_mask_clust(mask);


[mask , ii]=w_mask_fillin_once(mask , 1);
jj=ii;
if ii>0
    [mask , ii]=w_mask_fillin_once(mask ,1);
    jj=jj+ii;
    if ii>0
        [mask , ii]=w_mask_fillin_once(mask , 1);
        jj=jj+ii;
    end
end

mask = imfill(mask,'holes');

function [mask , nfill]= w_mask_fillin_once( mask , nside)
[nx , ny , nz]=size(mask);

nsx = nside ; nsy = nside ; nsz = nside;

mask=reshape(mask ,[] , 1);
nfill=0;

nxy=nx*ny;

nx2  =2*nx  ; nx3  = 3*nx  ; nx4  = 4*nx;
nxy2 =2*nxy ; nxy3 = 3*nxy ; nxy4 = 4*nxy;


for kk=nsz:(nz-nsz-1)
    kv=kk * nxy;
    for jj=nsy:(ny-nsy-1)
        jv=jj * nx + kv;
        for ii=nsx:(nx-nsx-1)
            iv = ii + jv;
            if ~mask(iv)
                %X-plant
                tag='';
                switch(nsx)
                    case 1
                        if mask(iv+1) && mask(iv-1)
                            mask(iv)=1;
                            tag='NextVox';
                            nfill=nfill+1;
                        end
                    case 2
                        if ((mask(iv+1) || mask(iv+2)) &&...
                            mask(iv-2) || mask(iv-2))
                           mask(iv)=1;
                           tag='NextVox';
                           nfill=nfill+1;
                        end
                    case 3
                        if (mask(iv+1) || mask(iv+2) ||...
                                mask(iv+3) &&...
                                mask(iv-1) || mask(iv-2) ||...
                                mask(iv-3))
                            mask(iv)=1;
                            tag='NextVox';
                            nfill=nfill+1;
                        end
                    case 4
                        if (mask(iv+1) || mask(iv+2) ||...
                            mask(iv+3) || mask(iv+4) &&...
                            mask(iv-1) || mask(iv-2) ||...
                            mask(iv-3) || mask(iv-4) )
                            mask(iv)=1;
                            tag='NextVox';
                            nfill=nfill+1;
                        end
                    otherwise
                end
                %Y-plant
                if ~strcmp(tag , 'NextVox')
                    switch(nsy)
                        case 1
                            if (mask(iv+nx) && mask(iv-nx))
                               mask(iv)=1;
                               tag='NextVox';
                               nfill=nfill+1;
                            end
                        case 2
                            if (mask(iv+nx) || mask(iv+nx2)  &&...
                                    mask(iv-nx) || mask(iv-nx2))
                                mask(iv)=1;
                                tag='NextVox';
                                nfill=nfill+1;
                            end
                        case 3
                            if (mask(iv+nx) || mask(iv+nx2) ||...
                                    mask(iv+nx3) && ...
                                    mask(iv-nx) || mask(iv-nx2) || ...
                                    mask(iv-nx3))
                                mask(iv)=1;
                                tag='NextVox';
                            end
                        case 4
                            if (mask(iv+nx) || mask(iv+nx2) || ...
                                    mask(iv+nx3) || mask(iv+nx4) &&...
                                    mask(iv-nx)  || mask(iv-nx2) ||...
                                    mask(iv-nx3) || mask(iv-nx4))
                                mask(iv)=1;
                                tag='NextVox';
                                nfill=nfill+1;
                            end
                    end
                end
                %Z-plant
                if ~strcmp(tag , 'NextVox')
                    switch(nsz)
                        case 1
                            if (mask(iv+nxy) && mask(iv-nxy))
                               mask(iv)=1;
                               
                               nfill=nfill+1;
                            end
                        case 2
                            if (mask(iv+nxy) || mask(iv+nxy2)  &&...
                                    mask(iv-nxy) || mask(iv-nxy2))
                                mask(iv)=1;
                                
                                nfill=nfill+1;
                            end
                        case 3
                            if (mask(iv+nxy) || mask(iv+nxy2) ||...
                                    mask(iv+nxy3) && ...
                                    mask(iv-nxy) || mask(iv-nxy2) || ...
                                    mask(iv-nxy3))
                                mask(iv)=1;
                                
                                nfill=nfill+1;
                            end
                        case 4
                            if (mask(iv+nxy) || mask(iv+nxy2) || ...
                                    mask(iv+nxy3) || mask(iv+nxy4) &&...
                                    mask(iv-nxy)  || mask(iv-nxy2) ||...
                                    mask(iv-nxy3) || mask(iv-nxy4))
                                mask(iv)=1;
                                
                                nfill=nfill+1;
                            end
                    end
                end
            end
        end
    end
end

mask=reshape(mask , [nx , ny , nz]);

function [mask , nfill]=w_mask_fillin_completely( mask , nside )

nfill=0;
[mask , kfill]=w_mask_fillin_once(mask , nside);
while( kfill>0 )
    [mask , kfill]=w_mask_fillin_once( mask , nside);
    nfill = nfill+kfill;
end
        
function mask =w_mask_clust(mask, RMM)
if ~exist('RMM', 'var')
    RMM=18;
end

mask=logical(mask);

[L , num]=bwlabeln(mask , RMM);
if num>=1
    MaxClusterSize=0;
    for i=1:num
        theCurrentCluster = L==i;
        MaxClusterSize = ...
            max([MaxClusterSize,length(find(theCurrentCluster))]);
    end
    
    for i=1:num
        theCurrentCluster = L==i;
        if MaxClusterSize==length(find(theCurrentCluster))
            mask_L=i;
        end
    end    
    mask = L==mask_L;
end

function mask=w_mask_erodemany( mask , npeel , peelthr)
[nx , ny , nz]=size(mask);
mask = reshape(mask , [] , 1);
nnn=zeros(size(mask));

nxy = nx*ny ; nxyz = nxy * nz;
for pp=1:npeel
    bpp = pp;
    for kk=0:(nz-1)
        kz = kk * nxy ; km = kz - nxy ; kp = kz + nxy;
        if kk == 0
            km = kz;
        end
        
        if kk == nz - 1;
            kp = kz;
        end
        
        for jj=0:(ny-1)
            jy = jj * nx ; jm = jy - nx ; jp = jy + nx;
            if jj == 0
                jm = jy;
            end
            
            if jj == ny-1
                jp = jy;
            end
            
            for ii=0:(nx-1)
                if mask( ii + jy + kz + 1)
                    im = ii - 1;
                    ip = ii + 1;
                    if ii == 0
                        im = 0;
                    end
                    
                    if ii == nx-1
                        ip = ii;
                    end
                
                    num = mask(im + jy + km + 1)...
                        + mask(ii + jm + km + 1) + mask(ii + jy + km + 1) + mask(ii + jp + km + 1)...
                        ...%    
                        + mask(ip + jy + km + 1)...
                        + mask(im + jm + kz + 1) + mask(im + jy + kz + 1) + mask(im + jp + kz + 1)...
                        ...%
                        + mask(ii + jm + kz + 1)                          + mask(ii + jp + kz + 1)...
                        + mask(ip + jm + kz + 1) + mask(ip + jy + kz + 1) + mask(ip + jp + kz + 1)...
                        ...%
                        + mask(im + jy + kp + 1)...
                        + mask(ii + jm + kp + 1) + mask(ii + jy + kp + 1) + mask(ii + jp + kp + 1)...
                        + mask(ip + jy + kp + 1) ;
                    
                    if num < peelthr 
                        nnn(ii + jy + kz + 1) = bpp;
                    end
                end
            end
        end
    end
    
    mask=mask.*(nnn==0);
end
qqq=zeros(size(mask));
for pp=npeel:-1:1
    bpp = pp;
    bth = ~(pp == npeel);
    
    for kk=0:(nz-1)
        kz = kk * nxy ; km = kz - nxy ; kp = kz + nxy ; 
        if kk==0
            km = kz ;
        end
        
        if kk==nz-1
            kp = kz ;
        end
        
        for jj=0:(ny-1)
            jy = jj*nx ; jm = jy - nx ; jp = jy+nx;
            if jj==0
                jm = jy ;
            end
            
            if jj==ny-1 
                jp = jy ;
            end
            
            for ii=0:(nx-1)
                if (nnn(ii + jy + kz + 1) >= bpp && ...
                        ~mask(ii + jy + kz + 1))
                    im = ii - 1;
                    ip = ii + 1;
                    if ii==0
                        im = 0;
                    end
                    
                    if ii==nx-1;
                        ip = ii;
                    end
                    qqq(ii + jy + + kz + 1) = mask(im + jy + + km + 1)...
                        + mask(ii + jm + km + 1) + mask(ii + jy + km + 1)...
                        + mask(ii + jp + km + 1) + mask(ip + jy + km + 1)...
                        ...%
                        + mask(im + jm + kz + 1) + mask(im + jy + kz + 1)... 
                        + mask(im + jp + kz + 1) + mask(ii + jm + kz + 1)...
                        ...%
                        + mask(ii + jp + kz + 1) + mask(ip + jm + kz + 1)...
                        + mask(ip + jy + kz + 1) + mask(ip + jp + kz + 1)...
                        ...%
                        + mask(im + jy + kp + 1) + mask(ii + jm + kp + 1)...
                        + mask(ii + jy + kp + 1) + mask(ii + jp + kp + 1)...
                        + mask(ip + jy + kp + 1) ;
                end
            end
        end
    end
    
    for ii=1:nxyz
        if qqq(ii) > bth
            mask(ii) = 1;
        end
    end
end

mask = reshape(mask , [nx , ny , nz]);

function clip_mask=CliplevelGradual(mfrac , mask)
[nx , ny , nz]=size(mask);

%Get the Center of mass of a 3D image (in index coordinates)
%Like AFNI's mri_get_cmass_3D
nxy = nx * ny;

mask1D=reshape(mask , [] , 1);
xx = 0 ; yy = 0 ; zz = 0 ;
sum = 0 ;
for kk=0:(nz-1)
    koff = kk * nxy;
    for jj=0:(ny-1)
        joff = koff + jj*nx;
        for ii=0:(nx-1)
            val = abs(mask1D(ii+joff+1));
            sum = sum + val;
            xx  = xx  + val*ii;
            yy  = yy  + val*jj;
            zz  = zz  + val*kk;
        end
    end
end

if sum > 0
    xx = xx/sum;
    yy = yy/sum;
    zz = zz/sum;
else
    xx = 0.5*(nx-1);
    yy = 0.5*(ny-1);
    zz = 0.5*(nz-1);
end

ic = round(xx);
jc = round(yy);
kc = round(zz);

clip=w_ClipLevel( mfrac , mask);
val = 0.333 * clip;

ii = round(0.01 * nx);
if ii<1
    ii=1;
end

jj = round(0.01 * ny); 
if jj<1
    jj=1;
end

kk = round(0.01 * nz); 
if kk<1
    kk=1;
end

icm = ic - ii + 1; icp = ic + ii + 1; 
if icm < 1
    icm = 1;
end

if icp > nx;
    icp = nx;
end

jcm = jc - jj + 1; jcp = jc + jj + 1; 
if jcm < 1
    jcm = 1;
end

if jcp > ny;
    jcp = ny;
end

kcm = kc - kk + 1; kcp = kc + kk + 1; 
if kcm < 1
    kcm = 1;
end

if kcp > nz;
    kcp = nz;
end

clip_000=w_cliplevel_partial(mask , mfrac , ...
    1   , icp , ...
    1   , jcp , ...
    1   , kcp);

clip_100=w_cliplevel_partial(mask , mfrac , ...
    icm , nx , ...
    1   , jcp , ...
    1   , kcp);

clip_010=w_cliplevel_partial(mask , mfrac , ...
    1   , icp , ...
    jcm , ny , ...
    1   , kcp);

clip_110=w_cliplevel_partial(mask , mfrac , ...
    icm , nx , ...
    jcm , ny , ...
    1   , kcp);

clip_001=w_cliplevel_partial(mask , mfrac , ...
    1   , icp , ...
    1   , jcp , ...
    kcm , nz);

clip_101=w_cliplevel_partial(mask , mfrac , ...
    icm , nx , ...
    1   ,  jcp, ...
    kcm , nz);

clip_011=w_cliplevel_partial(mask , mfrac , ...
    1   , icp , ...
    jcm , ny , ...
    kcm , nz);

clip_111=w_cliplevel_partial(mask , mfrac , ...
    icm , nx , ...
    jcm , ny , ...
    kcm , nz);
%1
if clip_000 < val
    clip_000 = val;
end
cv.clip_000=clip_000;
clip=clip_000;
clear clip_000
%2
if clip_100 < val
    clip_100 = val;
end
cv.clip_100=clip_100;
clip=[clip;clip_100];
clear clip_100
%3
if clip_010 < val
    clip_010 = val;
end
cv.clip_010=clip_010;
clip=[clip;clip_010];
clear clip_010
%4
if clip_110 < val
    clip_110 = val;
end
cv.clip_110=clip_110;
clip=[clip;clip_110];
clear clip_110
%5
if clip_001 < val
    clip_001 = val;
end
cv.clip_001=clip_001;
clip=[clip;clip_001];
clear clip_001
%6
if clip_101 < val
    clip_101 = val;
end
cv.clip_101=clip_101;
clip=[clip;clip_101];
clear clip_101
%7
if clip_011 < val
    clip_011 = val;
end
cv.clip_011=clip_011;
clip=[clip;clip_011];
clear clip_011
%8
if clip_111 < val
    clip_111 = val;
end
cv.clip_111=clip_111;
clip=[clip;clip_111];
clear clip_111

cv.clip_min=min(clip);
cv.clip_max=max(clip);

cv.x0 = 0.5 * ic ; cv.x1 = 0.5*(ic+nx-1) ;
cv.y0 = 0.5 * jc ; cv.y1 = 0.5*(jc+ny-1) ;
cv.z0 = 0.5 * kc ; cv.z1 = 0.5*(kc+nz-1) ;

if cv.x1 > cv.x0
    cv.dxi = 1/(cv.x1-cv.x0);
else
    cv.dxi = 0;
end

if cv.y1 > cv.y0
    cv.dyi = 1/(cv.y1-cv.y0);
else
    cv.dyi = 0;
end

if cv.z1 > cv.z0
    cv.dzi = 1/(cv.z1-cv.z0);
else
    cv.dzi = 0;
end

clip_mask=zeros(size(mask1D));

ijk=1;
for kk=0:nz-1
    for jj=0:ny-1
        for ii=0:nx-1
            clip_mask(ijk) = pointclip( ii , jj , kk ,cv);
            ijk = ijk + 1;
        end
    end
end

clip_mask=reshape(clip_mask , [nx , ny , nz]);
            

function clip=w_cliplevel_partial( mask , mfrac ,...
    xa , xb ,...
    ya , yb ,...
    za , zb)

tmp_mask = mask(xa:xb , ya:yb , za:zb);
clip=w_ClipLevel(mfrac , tmp_mask);
    
function clip=pointclip( ii , jj , kk , cv)

x1 = ( ii-cv.x0 ) * cv.dxi;
if x1 < 0
    x1 = 0;
elseif x1 > 1
    x1 = 1;
end

y1 = ( jj-cv.y0 ) * cv.dyi;
if y1 < 0
    y1 = 0;
elseif y1 > 1
    y1 = 1;
end

z1 = ( kk-cv.z0 ) * cv.dzi;
if z1 < 0
    z1 = 0;
elseif z1 > 1
    z1 = 1;
end

x0 = 1 - x1 ; y0 = 1 - y1 ; z0 = 1 - z1 ;

clip= cv.clip_000 * x0*y0*z0 + cv.clip_100 * x1*y0*z0...
    + cv.clip_010 * x0*y1*z0 + cv.clip_110 * x1*y1*z0...
    + cv.clip_001 * x0*y0*z1 + cv.clip_101 * x1*y0*z1...
    + cv.clip_011 * x0*y1*z1 + cv.clip_111 * x1*y1*z1;