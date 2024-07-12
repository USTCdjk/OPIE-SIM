clear all
close all
clc
addpath('.\OtfProvider\');
addpath('.\backgroundsuppression\');
addpath('.\readdata\');
addpath('.\SIMfunction\');
addpath('.\showpicture\');
addpath('.\fourierShift\');
addpath('.\WienerFilter\');
addpath('.\savefunction\');
addpath('.\lib\bfmatlab\');
addpath('.\lib\diplib\share\DIPimage');
%% Loading system OTF file
rawpixelsize = [65 65 100]; % pixel size and focal stack spacing (nm) 3D:85
NA = 1.49;                   % objective lens NA 1.49 3D:1.49  
refmed = 1.47;              % refractive index medium
refcov = 1.512;             % refractive index cover slip
refimm = 1.515;             % refractive index immersion medium
exwavelength =405;         % excitation wavelengths 488 561 638 668
emwavelength =525;         % emission wavelengths   525 610 665 680
fwd = 140e3;                % free working distance from objective to cover slip
depth = 0;                  % depth of imaged slice from the cover slip
numsteps=3;                 %样本相位个数
nrBands=2;                  %相干光数量，例如2DSIM，S（k）与S（k±p）band就是两个，3D，三个
numangles=3;                %结构光方向数
readmode=2;%1.paz 2.pza(OMX) 

% readplane=[1,2,3,4,5,6,7,8,9,10];
% readplane=[1,5,9];
readplane=[1];
ok=0;
use_opt=1;
freq_initial=[];
usebackgroundsupression=0;
show_pic=0;
if(ok==1)
    freq_initial=[];
end
%% 3相位2D样本，可通过样本中同numangles，Z的图片个数来判断
% filepath='My_testpic_3P\RawSIMData_gt-2.tif';
% filepath='My_testpic_3P\RawSIMData2_level_05.tif';
% filepath='My_testpic_3P\RawSIMData_gt-1.tif';
% filepath='My_testpic_3P\638\multicolordata-638.tif';
% filepath='My_testpic_3P\488\Stack488P3.tif';
filepath='论文\效果2\405\max\活细胞多色405.tif';
% filepath='论文\效果2\405\t=6\405.tif';
% filepath='My_testpic_3P\Raw_BF_488_actin.tif';
%% 5相位2D样本
% filepath='My_testpic_5P\OMX_LSEC_Actin_525nm.tif';
% filepath='My_testpic_5P\OMX_LSEC_Membrane_680nm.tif';
% filepath='My_testpic_5P\OMX_Tetraspeck200_680nm.tif';
% filepath='My_testpic_5P\488\488stack.tif';
% filepath='My_testpic_5P\638\t8_p5_638_z100_data1.tif';
% filepath='论文\simulate_3DMicrotubules.tif';
% filepath='My_testpic_5P\405\405_t4.6_z100.tif';
% [data2D,expand_value]=read3Ddata2(filepath,numsteps,numangles);  %%(Nx,Ny,numsteps,Nz,angle)  %读取图片到数组
[data2D,expand_value,focu_plane]=readdata(filepath,numsteps,numangles,readmode); %%(Nx,Ny,numsteps,Nz,angle)  %读取图片到数组
close all;
% if(usebackgroundsupression==1)
%     [data2D]=backgroundsuppression(data2D,exwavelength,[1,50]);%%第三个矩阵代表BF强度10/50值越小BF强度越大
% end
%% 3相位裁剪后的2D样本
% filepath='testpic_unzip\snap_01_000';
% filepath='testpic_unzip\Actin_1-1_z001_c00';
% [data2D,expand_value]=read2Ddata_unzip(filepath,numsteps,numangles);%(Nx,Ny,numsteps,Nz,angle)Nx为行数

%% 读取尼康样本
% filepath='nikang\sim01z4.tif';
% [data2D,expand_value]=read2Ddata_nikang(filepath,numsteps,numangles,512,512,1);%(Nx,Ny,numsteps,Nz,angle)  %读取图片到数组

rawpixelsize(1) = rawpixelsize(1)*expand_value.x;
rawpixelsize(2) = rawpixelsize(2)*expand_value.y;

if(show_pic==1)
    showall(data2D);%按照相位，层数，角度显示所有图片
    close all;
end

notchwidthxy1 = 0.4;    % notch_width to design Filter
notchdips1 = 0.98;              % notch_depth to design Filter
notchwidthxy2 = 0.4;    % notch_width to notch
notchdips2 = 0.98;              % notch_depth to notch
OTFflag = 1;                    % if OTFflag == 1,simulated OTF; if OTFfla == 0 ,experimental OTF
OTF_name = '.\OMX-OTFs\OMX-OTF-528nm-2d.xml';
attenuation = 1;                %阻尼因子

[numpixelsx,numpixelsy,numsteps,numfocus,numangles] = size(data2D.allimages_in_Zstack);


for fl=1:length(readplane)
        Z_now=readplane(fl);
        if(Z_now<numfocus+1 && Z_now>0)
            data2D.allimages_in= data2D.allimages_in_Zstack(:,:,:,Z_now,:);
            [data2D,freq,ang,pha,module0,final_image,fair_image]=My_process_data(data2D,nrBands,rawpixelsize,NA,refmed,refcov,refimm,exwavelength,...
                emwavelength,fwd,depth,notchwidthxy1,notchdips1,notchwidthxy2,notchdips2,ok,OTFflag,OTF_name,attenuation,freq_initial,usebackgroundsupression,use_opt);
%             ok=1;
            data2D.final_image(:,:,fl)=final_image;
            disp('Z_now='+string(Z_now))
            close all
        end
end
[final_pic]=generateresult(data2D,usebackgroundsupression);

clear all
close all
% save('.\save\data2D.mat','data2D');