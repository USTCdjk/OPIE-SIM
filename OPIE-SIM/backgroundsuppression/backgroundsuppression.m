function [data2D]=backgroundsuppression(data2D,exwavelength,depth)
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
[imgx_psf,imgy_psf,imgz_psf]=size(PSF_3D);
PSF_3D=PSF_3D./sum(sum(PSF_3D(:)));
[imgx1,imgy1,numsteps,Nz,numangles]=size(data2D.allimages_in_Zstack);
if((imgx1<imgx_psf)||(imgy1<imgy_psf))
    PSF_3D=resizepsf(PSF_3D,[imgx1,imgy1]);
end
OTF_z=zeros(size(PSF_3D));
for ii=1:imgz_psf
    OTF_z(:,:,ii)=abs(Ft2(PSF_3D(:,:,ii)));
end
OTF_sum_infocus=sum(OTF_z(:,:,floor(imgz_psf/2)+1-depth_in:floor(imgz_psf/2)+1+depth_in),3);
OTF_sum_outfocus=sum(OTF_z(:,:,floor(imgz_psf/2)-depth_out+1:floor(imgz_psf/2)-depth_in),3)+sum(OTF_z(:,:,floor(imgz_psf/2)+depth_in+2:floor(imgz_psf/2)+depth_out+1),3);

err=0.12;
OTF_l_mask=OTF_sum_outfocus./((OTF_sum_infocus+OTF_sum_outfocus)+err);
OTF_l_mask=resizemask(OTF_l_mask,[imgx1,imgy1]);

blur_img_defocus_new=zeros(imgx1,imgx1,Nz);

for angle=1:numangles
    for j=1:Nz
        for phase=1:numsteps
            bg=abs(iFt2(Ft2(data2D.allimages_in_Zstack(:,:,phase,j,angle)).*OTF_l_mask));
            blur_img_defocus_new(:,:,phase+numsteps*(j-1)+numsteps*Nz*(angle-1))=data2D.allimages_in_Zstack(:,:,phase,j,angle)-bg;
            data2D.allimages_in_Zstack(:,:,phase,j,angle)=blur_img_defocus_new(:,:,phase+numsteps*(j-1)+numsteps*Nz*(angle-1));
        end
    end
end

imwritestack_16(blur_img_defocus_new,['save\','Raw_BF_image.tif'])

function [out] = resizepsf( in,tosize )
    siz=size(in);
    w=siz(2);
    h=siz(1);
    out=in(h/2-tosize(1)/2+1:tosize(1)/2+h/2,w/2-tosize(2)/2+1:tosize(2)/2+w/2,:);
end

function [out] = resizemask( in,tosize )
    siz=size(in);
    w=siz(2);
    h=siz(1);
    out=zeros(tosize(1),tosize(2));
    out(tosize(1)/2-h/2+1:tosize(1)/2+h/2,tosize(2)/2-w/2+1:tosize(2)/2+w/2)=in;
end

end