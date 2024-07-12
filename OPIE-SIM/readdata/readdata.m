function [data2D,expand_value,focu_plane]=readdata(filepath,numsteps,numangles,readmode,cutsize,offset)
data2D.datasets = filepath;
data2D.numsteps = numsteps;
data2D.numangles = numangles;


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
if(nargin<5)
   picsize_x=Xi;
   picsize_y=Yi;
else
    picsize_x=cutsize(1);
    picsize_y=cutsize(2);
end
if(nargin<6)
  offset_x=0;
  offset_y=0;
else
    offset_x=offset(1);
    offset_y=offset(2);
end

cutpixelsx=[1+offset_x,picsize_x+offset_x];%行
cutpixelsy=[1+offset_y,picsize_y+offset_y];%列
% cutpixelsx=[100,200];
% cutpixelsy=[100,200];
data2D.numpixelsx=picsize_x;
data2D.numpixelsy=picsize_y;
data2D.imgSize=data2D.numpixelsx;

%截图后放大相当于放大了实域坐标，则频域坐标被缩小
Image=zeros(Yi,Xi,SL);
Nz=SL/(numsteps*numangles);
focu_plane=ceil(Nz/2);

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
%     ishow(Image(:,:,i),'');
end

if (cutpixelsx(2)-cutpixelsx(1)+1)~=Xi ||(cutpixelsy(2)-cutpixelsy(1)+1)~=Yi
    expand_value.x=(cutpixelsx(2)-cutpixelsx(1)+1)/data2D.numpixelsx;
    expand_value.y=(cutpixelsy(2)-cutpixelsy(1)+1)/data2D.numpixelsy;
else
    expand_value.x=Xi/data2D.numpixelsx;
    expand_value.y=Yi/data2D.numpixelsy;
end

if(readmode==1)
    for j=1:Nz
        for angle=1:numangles
            for phase=1:numsteps
                if (cutpixelsx(2)-cutpixelsx(1)+1)~=Xi ||(cutpixelsy(2)-cutpixelsy(1)+1)~=Yi
                    a = imresize(double( Image(cutpixelsx(1):cutpixelsx(2),cutpixelsy(1):cutpixelsy(2),phase+numsteps*(angle-1)+numsteps*numangles*(j-1))),[data2D.numpixelsx,data2D.numpixelsy],'bilinear');
                else
                    a = imresize(double( Image(1:end,1:end,phase+numsteps*(angle-1)+numsteps*numangles*(j-1))),[data2D.numpixelsx,data2D.numpixelsy],'bilinear');
                    a = imresize(a,[data2D.numpixelsx,data2D.numpixelsy],'bilinear');
                end
                data2D.allimages_in_Zstack(:,:,phase,j,angle)= a;
            end
        end
    end
elseif(readmode==2)
    for angle=1:numangles
        for j=1:Nz
            for phase=1:numsteps
                if (cutpixelsx(2)-cutpixelsx(1)+1)~=Xi ||(cutpixelsy(2)-cutpixelsy(1)+1)~=Yi
                    a = imresize(double( Image(cutpixelsx(1):cutpixelsx(2),cutpixelsy(1):cutpixelsy(2),phase+numsteps*(j-1)+numsteps*Nz*(angle-1))),[data2D.numpixelsx,data2D.numpixelsy],'bilinear');
                else
                    a = imresize(double( Image(1:end,1:end,phase+numsteps*(j-1)+numsteps*Nz*(angle-1))),[data2D.numpixelsx,data2D.numpixelsy],'bilinear');
                    a = imresize(a,[data2D.numpixelsx,data2D.numpixelsy],'bilinear');
                end
%                 disp(phase+numsteps*(j-1)+numsteps*Nz*(angle-1));
                data2D.allimages_in_Zstack(:,:,phase,j,angle)= a;
%                 ishow(a,'');
            end
        end
    end
end
