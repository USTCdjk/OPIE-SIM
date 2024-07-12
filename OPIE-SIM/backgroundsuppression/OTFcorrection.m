function [outimg]=OTFcorrection(fftimg,exwavelength,numsteps,numangles,depth)
%% background suppression
Ft2 = @(x) fftshift(fft2(ifftshift(x)));
iFt2 = @(x) fftshift(ifft2(ifftshift(x)));
if(nargin<3)
    depth_in=4;
    depth_out=50;
else
    depth_in=depth(1);
    depth_out=depth(2);
end
PSF_3D=imreadstack('PSF RW_lambda='+string(exwavelength)+'.tif');
% PSF_3D=imreadstack('PSF BW.tif');
[imgx_psf,imgy_psf,imgz_psf]=size(PSF_3D);
PSF_3D=PSF_3D./sum(sum(PSF_3D(:)));
[imgx1,imgy1,~]=size(fftimg);
if((imgx1<imgx_psf)||(imgy1<imgy_psf))
    PSF_3D=resizepsf(PSF_3D,[imgx1,imgy1]);
end
OTF_z=zeros(size(PSF_3D));
for ii=1:imgz_psf
    OTF_z(:,:,ii)=abs(Ft2(PSF_3D(:,:,ii)));
end
OTF_sum_infocus=sum(OTF_z(:,:,floor(imgz_psf/2)+1-depth_in:floor(imgz_psf/2)+1+depth_in),3);
OTF_sum_outfocus=sum(OTF_z(:,:,floor(imgz_psf/2)-depth_out+1:floor(imgz_psf/2)-depth_in),3)+sum(OTF_z(:,:,floor(imgz_psf/2)+depth_in+2:floor(imgz_psf/2)+depth_out+1),3);
w1=0.1; %0.1
w2=0.8; %0.8
OTF_l_mask=(OTF_sum_infocus+w1)./((OTF_sum_infocus+OTF_sum_outfocus)+w2);
if((imgx1>imgx_psf)||(imgy1>imgy_psf))
OTF_l_mask=resizemask(OTF_l_mask,[imgx1,imgy1]);
end
outimg=zeros(imgx1,imgx1,numangles*numsteps);
blur_img_defocus_new=zeros(imgx1,imgx1,numangles*numsteps);
for i=1:numangles*numsteps
    bg=abs(iFt2(Ft2(fftimg(:,:,i)).*OTF_l_mask));
    blur_img_defocus_new(:,:,i)=bg;
    outimg(:,:,i)=blur_img_defocus_new(:,:,i);
end


function [out] = resizepsf( in,tosize )
    siz=size(in);
    w=siz(2);
    h=siz(1);
    out=in(h/2-tosize(1)/2+1:tosize(1)/2+h/2,w/2-tosize(2)/2+1:tosize(2)/2+w/2,:);
end

function [out] = resizemask( in,tosize )
    out=imresize(double(in),tosize,'bilinear');
end

end