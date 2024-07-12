function [data2D] = save_SIM(FsumA,FperiA,Fcent,data2D,Z_now)
%% recontructed SIM images
data2D.WF_FR(:,:,Z_now)=Fcent;
data2D.offcenter_FR(:,:,Z_now)=FperiA;
data2D.SIM_FR(:,:,Z_now)=FsumA;

% Dsum = real( ifft2(fftshift(Fsum)) );
% Dperi = real( ifft2(fftshift(Fperi)) );
% %% 最后解卷积
% data2D.SIM_AfterDe(:,:,Z_now)=deconvlucy(Dsum,data2D.PSF_shiftfinal_SIM,5);
% data2D.offcenter_AfterDe(:,:,Z_now)=deconvlucy(Dperi,data2D.PSF_shiftfinal_offcenter,5);