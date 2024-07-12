function [data3D,expand_value]=read3Ddata2(filepath,numsteps,numangles)
data3D.datasets = filepath;
data3D.numsteps = numsteps;
data3D.numangles = numangles;


%% read dv/tiff/tif
info=imfinfo(filepath);
tif='tif';
format=info.Format;
if(strcmp(format,tif)==0)
    disp('error');
end
%Wi代表矩阵dim=1，He代表矩阵dim=2
SL=size(info,1);
Yi=info.Width;%图像的宽度在matlab里dim=2
Xi=info.Height;%图像的高度在matlab里dim=1
offset_x=0;
offset_y=0;
cutpixelsx=[1+offset_x,1024+offset_x];%行
cutpixelsy=[1+offset_y,1024+offset_y];%列
% cutpixelsx=[100,200];
% cutpixelsy=[100,200];
data3D.numpixelsx=1024;
data3D.numpixelsy=1024;
data3D.imgSize=data3D.numpixelsx;

%截图后放大相当于放大了实域坐标，则频域坐标被缩小
Image=zeros(Yi,Xi,SL);
Nz=SL/(numsteps*numangles);
data3D.Nz=Nz;

if (cutpixelsx(2)>Xi)
    cutpixelsx(2)=Xi;
end
if (cutpixelsx(1)<1) 
    cutpixelsx(1)=1;
    cutpixelsy(1)=1;
end
if (cutpixelsy(2)>Yi)
    cutpixelsy(2)=Yi;
end
if (cutpixelsy(1)<1)
    cutpixelsy(1)=1;
end


for i=1:SL
    Image(:,:,i)=imread(filepath,i);
end

if (cutpixelsx(2)-cutpixelsx(1)+1)~=Xi ||(cutpixelsy(2)-cutpixelsy(1)+1)~=Yi
    expand_value.x=(cutpixelsx(2)-cutpixelsx(1)+1)/data3D.numpixelsx;
    expand_value.y=(cutpixelsy(2)-cutpixelsy(1)+1)/data3D.numpixelsy;
else
    expand_value.x=Xi/data3D.numpixelsx;
    expand_value.y=Yi/data3D.numpixelsy;
end


for j=1:Nz
    for angle=1:numangles
        for phase=1:numsteps
            if (cutpixelsx(2)-cutpixelsx(1)+1)~=Xi ||(cutpixelsy(2)-cutpixelsy(1)+1)~=Yi
                a = imresize(double( Image(cutpixelsx(1):cutpixelsx(2),cutpixelsy(1):cutpixelsy(2),phase+numsteps*(j-1)+numsteps*Nz*(angle-1))),[data3D.numpixelsx,data3D.numpixelsy],'bilinear');
            else
                a = imresize(double( Image(1:end,1:end,phase+numsteps*(angle-1)+numsteps*numangles*(j-1))),[data3D.numpixelsx,data3D.numpixelsy],'bilinear');
                a = imresize(a,[data3D.numpixelsx,data3D.numpixelsy],'bilinear');
            end
            data3D.allimages_in_Zstack(:,:,phase,j,angle)= a;
        end
    end
end

