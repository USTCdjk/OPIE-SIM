function [data2D,freq,ang,pha,module0,final_image,fair_image] = My_process_data(data2D,nrBands,rawpixelsize,NA,refmed,refcov,refimm,exwavelength,...
    emwavelength,fwd,depth,notchwidthxy1,notchdips1,notchwidthxy2,notchdips2,ok,OTFflag,OTF_name,attenuation,freq_initial,usebackgroundsupression,use_opt)

%% Get some parameter
[numpixelsx,numpixelsy,numsteps,numfocus,numangles] = size(data2D.allimages_in);
[~,~,~,numfocus_sum,~] = size(data2D.allimages_in_Zstack);
data2D.rawpixelsize = rawpixelsize;
data2D.NA = NA;
data2D.refmed = refmed;
data2D.refcov = refcov;
data2D.refimm = refimm;  %浸润层反射率
data2D.exwavelength = exwavelength;
data2D.emwavelength = emwavelength;
data2D.fwd = fwd;
data2D.depth = depth;
data2D.xemit = 0;        % x-position focal point
data2D.yemit = 0;        % y-position focal point
data2D.zemit = 0;        % z-position focal point
data2D.attStrength = 0.9;
if(usebackgroundsupression==0)
    err=-0.01;
else
    err=0;
end
use_predec=1;
use_afterdec=0;
data2D.nrDirs=numangles;
data2D.nrBands = nrBands;
data2D.maxorder = 3;

% Below are the parameter to design the notchfilter, we have selected a
% proper value which fits most samples. But for better reconstruction,
% you can change them.
%
data2D.notchwidthxy1 = notchwidthxy1;
data2D.notchdips1 = notchdips1;
data2D.notchwidthxy2 = notchwidthxy2;
data2D.notchdips2 = notchdips2;


for jangle = 1:numangles
    tempimage = squeeze(sum(data2D.allimages_in(:,:,:,:,jangle),4));%1:4:13
%     tempimage = squeeze(data2D.allimages_in(:,:,:,ceil(numfocus/2),jangle));   %squeeze维度中只有一层的维度去掉，[numpixelsx,numpixelsy,numsteps,numfocus,numangles]
%     tempimage = weightCNR(tempimage,MCNR);
    Snoisy(:,:,(jangle-1)*numsteps+1:(jangle-1)*numsteps+numsteps) = tempimage;%只取中间Z的numangles*numsteps张图片
end





%% Get the frequency, phase and modulation depth
if(data2D.nrBands==2 && use_opt==1)
    [data2D,SeparateII_cell,find_peak_ok] = find_illumination_pattern2(Snoisy,data2D,ok,OTFflag,OTF_name,attenuation,use_predec,usebackgroundsupression);%输入为调制度最高的Nz的累加，此处就一层
elseif(data2D.nrBands==2 && use_opt==0)
    [data2D,SeparateII_cell] = find_illumination_pattern(Snoisy,data2D,ok,OTFflag,OTF_name,attenuation,use_predec,usebackgroundsupression);
elseif(data2D.nrBands==3)
    [data2D,SeparateII_cell] = find_illumination_pattern3(Snoisy,data2D,use_predec,usebackgroundsupression);
end

if(use_opt==1 && find_peak_ok==0)
    [data2D,SeparateII_cell] = find_illumination_pattern(Snoisy,data2D,ok,OTFflag,OTF_name,attenuation,use_predec,usebackgroundsupression);
end
freq = data2D.allpatternpitch; %三个方向上的波矢量
ang = data2D.allpatternangle;  %三个方向角度
module0= data2D.allmodule;%三个方向的调制度
pha = data2D.allpatternphases; %三个方向相位，列代表不同方向


clear Snoisy tempimage

if((data2D.nrBands==2))
    allmodule = [0 0]; %调制度求均值
    for jangle = 1:numangles
        allmodule = allmodule + data2D.allmodule(:,jangle)'/numangles;
    end

    disp(['    parameter of dir1 is: freq=',num2str(sqrt(freq(1,1)^2+freq(2,1)^2)),', ang=',num2str(ang(1,1)),', pha=',num2str(pha(1,1)),', module1=',num2str(data2D.allmodule(2,1))]);
    disp(['    parameter of dir2 is: freq=',num2str(sqrt(freq(1,2)^2+freq(2,2)^2)),', ang=',num2str(ang(1,2)),', pha=',num2str(pha(1,2)),', module1=',num2str(data2D.allmodule(2,2))]);
    disp(['    parameter of dir3 is: freq=',num2str(sqrt(freq(1,3)^2+freq(2,3)^2)),', ang=',num2str(ang(1,3)),', pha=',num2str(pha(1,3)),', module1=',num2str(data2D.allmodule(2,3))]);

    data2D.allmoduleave = allmodule;
    clear allmodule
elseif(data2D.nrBands==3)
    allmodule = [0 0 0];
    for jangle = 1:numangles
        allmodule = allmodule + data2D.allmodule(:,jangle)'/numangles;
    end

    disp(['    parameter of dir1 is: freq=',num2str(sqrt(freq(1,1)^2+freq(2,1)^2)),', ang=',num2str(ang(1,1)),', pha=',num2str(pha(1,1)),', module1=',num2str(data2D.allmodule(2,1)),', module2=',num2str(data2D.allmodule(3,1))]);
    disp(['    parameter of dir2 is: freq=',num2str(sqrt(freq(1,2)^2+freq(2,2)^2)),', ang=',num2str(ang(1,2)),', pha=',num2str(pha(1,2)),', module1=',num2str(data2D.allmodule(2,2)),', module2=',num2str(data2D.allmodule(3,2))]);
    disp(['    parameter of dir3 is: freq=',num2str(sqrt(freq(1,3)^2+freq(2,3)^2)),', ang=',num2str(ang(1,3)),', pha=',num2str(pha(1,3)),', module1=',num2str(data2D.allmodule(2,3)),', module2=',num2str(data2D.allmodule(3,3))]);

end




fftDirectlyCombined=zeros(numpixelsx*2,numpixelsy*2);
fftWF=zeros(numpixelsx*2,numpixelsy*2);
fftoffcenter=zeros(numpixelsx*2,numpixelsy*2);
OTF_shifted_combined=zeros(numpixelsx*2,numpixelsy*2);
OTFo=zeros(numpixelsx,numpixelsy);
OTFo=applyOtf(OTFo,data2D.OtfProvider,1,0,0,1,1);
% save_as_mat(OTFo,'OTFatt');
for I=1:numangles
    separate=SeparateII_cell{1,I} ; 
    
    shifted=zeros(2*numpixelsx,2*numpixelsy,data2D.maxorder,numangles);
    OTF_shifted=zeros(2*numpixelsx,2*numpixelsy,data2D.maxorder,numangles);
%     ishow(separate(:,:,1).*OTFo,'');
    shifted(:,:,1,I)=placeFreq(separate(:,:,1).*OTFo);
    ishow(shifted(:,:,1,I),'O');
    for b=2:data2D.nrBands
        pos=b*2-2;
        neg=b*2-1;
        shifted(:,:,pos,I)=placeFreq(separate(:,:,pos).*OTFo);
        shifted(:,:,neg,I)=placeFreq(separate(:,:,neg).*OTFo);
        
        
%         ishow(shifted(:,:,pos),'shifted(:,:,pos)');
%         ishow(shifted(:,:,neg),'shifted(:,:,neg)');
    end
    OTF_shifted(:,:,1,I)=applyOtf(shifted(:,:,1,I),data2D.OtfProvider,1,0,0,0,1);
%     shifted(:,:,1,I)=shifted(:,:,1,I).*OTF_shifted(:,:,1,I);
    for b=2:data2D.nrBands
        pos=b*2-2;
        neg=b*2-1;
        OTF_shifted(:,:,pos,I)=applyOtf(shifted(:,:,pos,I),data2D.OtfProvider,b,-(b-1)*freq(1,I),-(b-1)*freq(2,I),0,1);
        OTF_shifted(:,:,neg,I)=applyOtf(shifted(:,:,neg,I),data2D.OtfProvider,b,(b-1)*freq(1,I),(b-1)*freq(2,I),0,1);
%         shifted(:,:,pos,I)=shifted(:,:,pos,I).*OTF_shifted(:,:,1,I);
%         shifted(:,:,neg,I)=shifted(:,:,neg,I).*OTF_shifted(:,:,1,I);
%         ishow(shifted(:,:,pos,I),'shifted(:,:,pos)');
%         ishow(shifted(:,:,neg,I),'shifted(:,:,neg)');
        shifted(:,:,pos,I)=NfourierShift(shifted(:,:,pos,I),(-(b-1)*freq(1,I)),(-(b-1)*freq(2,I)));
        shifted(:,:,neg,I)=NfourierShift(shifted(:,:,neg,I),((b-1)*freq(1,I)),((b-1)*freq(2,I)));
        ishow(shifted(:,:,pos,I),'shifted(:,:,pos)');
        ishow(shifted(:,:,neg,I),'shifted(:,:,neg)');
    end
    
    for J=1:data2D.nrBands*2-1
        if(J==1)
           fftDirectlyCombined=fftDirectlyCombined+shifted(:,:,J,I)/1;
           fftWF=fftWF+shifted(:,:,J,I);
        else
           fftDirectlyCombined=fftDirectlyCombined+shifted(:,:,J,I);
           fftoffcenter=fftoffcenter+shifted(:,:,J,I);
        end
        OTF_shifted_combined=OTF_shifted_combined+OTF_shifted(:,:,J,I);
    end
end

% ishow(OTF_shifted_combined,'OTF_shifted_combined');
% ishow(OTF_shifted(:,:,1,I),'OTF_aS');
% PSF_aS=abs(otf2psf((OTF_shifted_combined)));
% PSF_bS=abs(otf2psf((OTF_shifted(:,:,1,I))));
% PSF_aS_FWHM=PSF_FWHM(PSF_aS);
% PSF_bS_FWHM=PSF_FWHM(PSF_bS);
% ishow(PSF_bS,'PSF_bS');
% ishow(PSF_aS,'PSF_aS');

% ishow(fftDirectlyCombined,'fftDirectlyCombined');
% ishow_Re(fftWF,'fftWF');
WF_image=real( ifft2(fftshift(fftWF)) );
WF_image=My_normalize(WF_image);
% if(err~=0)
%     WF_image(WF_image<err*max(max(WF_image)))=0;%%
% end
% save_tif_32(WF_image,'.\save\WF_image.tif')

% ishow_Re(fftoffcenter,'fftoffcenter');
offcenter_image=real( ifft2(fftshift(fftoffcenter)) );
offcenter_image=My_normalize(offcenter_image);
if(err~=0)
    offcenter_image(offcenter_image<err*max(max(offcenter_image)))=0;%%
end
% save_tif_32(offcenter_image,'.\save\offcenter_image.tif')

% ishow_Re(fftDirectlyCombined,'fftDirectlyCombined');
Directly_image=real( ifft2(fftshift(fftDirectlyCombined)) );
Directly_image=My_normalize(Directly_image);
if(err~=0)
    Directly_image(Directly_image<err*max(max(Directly_image)))=0;%%
end
% save_tif_32(Directly_image,'.\save\Directly_image.tif')
% ishow(OTF_shifted_combined,'shifted(:,:,neg)');
%% 在这里加入DPR算法，原图是Directly_image，PSF_aS_FWHM是SIM之后的PSF半峰全宽，PSF_bS_FWHM是SIM之前的PSF半峰全宽
%算法输入为Directly_image，PSF_aS_FWHM，生成结果要与fair_image和Hifi_Initial_image比较





w1_fair=0.3;
w1_Hifi=1.2;  % Initial optimization Wiener constant：[0.9-2.5]
w2_Hifi=0.5;
%% 
data2D.cutoff=1000/(0.5*data2D.lambda/data2D.NA);                        
data2D.sampleLateral=ceil(data2D.cutoff/data2D.cyclesPerMicron)+1;  
K0 = [sqrt(freq(1,1)^2+freq(2,1)^2),sqrt(freq(1,2)^2+freq(2,2)^2),sqrt(freq(1,3)^2+freq(2,3)^2)];

K=max([ceil(K0)]);
if data2D.nrBands==2
    cutoff=floor(1*K)/data2D.sampleLateral+1.0;
elseif	data2D.nrBands==3
    cutoff=floor(2*K)/data2D.sampleLateral+1.0;
end
otfHiFi=zeros(2*numpixelsx,2*numpixelsy);
otfHiFi=writeApoVector(otfHiFi,data2D.OtfProvider,cutoff);      % Ideal OTF
% save_as_mat(otfHiFi,'Apo');
Mask=zeros(2*numpixelsx,2*numpixelsy);


%% HiFi-SIM
% Step 1
% if size(Iraw,3)==9
if data2D.nrBands==2
    wFilter1=WienerFilterW1_2D(data2D);
else
    wFilter1=WienerFilterW1_3D(data2D);
end
Mask(wFilter1.wDenom~=0)=1;
Wk1_Hifi=otfHiFi./(wFilter1.wDenom+w1_Hifi^2);

fftInitialHiFi=fftDirectlyCombined.*Wk1_Hifi.*Mask;
% ishow(fftInitialHiFi,'fftHiFi');
% ishow_Re(fftInitialHiFi,'InitialHiFi');
Initial_image=real( ifft2(fftshift(fftInitialHiFi)) );
Initial_image=My_normalize(Initial_image);
% if(err~=0)
% Initial_image(Initial_image<0)=err*max(max(Initial_image));%%
% end
% save_tif_32(Initial_image,'.\save\Hifi_Initial_image.tif')

Wk_fair=otfHiFi./(wFilter1.wDenom+w1_fair^2);
% save_as_mat(Wk_fair.*Mask,'Wfilter');
fftfair=fftDirectlyCombined.*Wk_fair.*Mask;
% ishow_Re(fftfair,'fair');
fair_image=real( ifft2(fftshift(fftfair)) );
% fair_image=My_normalize(fair_image);
if(err~=0)
    fair_image(fair_image<err*max(max(fair_image)))=0;%%
end
% save_tif_32(fair_image,'.\save\fair_image.tif')

% Step 2
if data2D.nrBands==2
    wFilter2=WienerFilterW2_2D(data2D);
else
    wFilter2=WienerFilterW2_3D(data2D);
end


ApoFWHM=0.5*(cutoff-1);
ApoFWHM=min(0.5,round(ApoFWHM*100)/100);

apo= apodize_gauss([2*numpixelsx,2*numpixelsy], struct('rad',ApoFWHM));  
Wk2_Hifi=apo./(wFilter2.wDenom+w2_Hifi^2);
W_filter=Wk1_Hifi.*Wk2_Hifi.*Mask;
fftHiFi=fftDirectlyCombined.*W_filter;
save_tif_32(abs(fftHiFi).^(1/10),'.\save\fftHiFi.tif')
final_image=real( ifft2(fftshift(fftHiFi)) );
% final_image=My_normalize(final_image);
PSF_final=abs(otf2psf((otfHiFi)));
if(use_afterdec==1)
    final_image=deconvlucy(final_image,PSF_final,5);
end
if(err~=0)
    final_image(final_image<err*max(max(final_image)))=0;%%
end
ishow_Re(fftHiFi,'HiFi');
% save_tif_32(final_image,'.\save\Hifi_final_image.tif')
ishow(Wk1_Hifi.*Mask,'Wk1');
ishow(Wk2_Hifi.*Mask,'Wk2');
ishow(W_filter,'W_filter');
end