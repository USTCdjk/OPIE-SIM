function [dataparams,SeparateII_cell] = find_illumination_pattern3(Snoisy,dataparams,use_predec,usebackgroundsupression)

%% Basic parameter
numangles = dataparams.numangles;
numsteps = dataparams.numsteps;
numpixelsx = dataparams.numpixelsx;
numpixelsy = dataparams.numpixelsy;
NA = dataparams.NA;
lambda = dataparams.emwavelength;
dataparams.lambda = lambda;
cyclesPerMicron = double(1./(dataparams.numpixelsx*dataparams.rawpixelsize(1)*0.001));
dataparams.cyclesPerMicron =cyclesPerMicron;
%% Calculate OTF

OtfProvider = SimOtfProvider(dataparams,NA,lambda,1);
PSF_cac=abs(otf2psf((OtfProvider.otf)));
dataparams.OtfProvider=OtfProvider;

if(use_predec==1)
    Temp=importImages(Snoisy);
    IIraw=deconvlucy(Temp,PSF_cac,5);
else
    IIraw=importImages(Snoisy);
end
if(usebackgroundsupression==1)
    IIraw=OTFcorrection(IIraw,lambda,numsteps,numangles,[5,50]);
end
IIrawFFT = zeros(numpixelsx,numpixelsy,numsteps*numangles);
for I=1:numangles*numsteps
    IIrawFFT(:,:,I)=FFT2D(IIraw(:,:,I),false);
%     ishow( IIrawFFT(:,:,I),'N'); %显示每个图片傅里叶变换结果
end

% for jangle = 1:numangles
%     for nstep = 1:numsteps
%         save_as_mat(IIrawFFT(:,:,(jangle-1)*numsteps+nstep),'pic_afterFFT_phase'+string(nstep)+'_dir'+string(jangle));
%     end
% end

cnt=[numpixelsx/2+1,numpixelsy/2+1];
dataparams.cutoff=1000/(0.5*dataparams.lambda/dataparams.NA);
[x,y]=meshgrid(1:numpixelsx,1:numpixelsy);
rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);
Mask=double(rad<=1.0*(dataparams.cutoff/dataparams.cyclesPerMicron+1));
NotchFilter0=getotfAtt(numpixelsx,cyclesPerMicron,0.5*dataparams.cutoff,0,0);
NotchFilter=NotchFilter0.*Mask;
% ishow(NotchFilter,'NotchFilter');

Mask2=double(rad<=1.10*(dataparams.cutoff/dataparams.cyclesPerMicron+1));
NotchFilter2=NotchFilter0.*Mask2;

CrossCorrelation=zeros(size(Mask2,1),size(Mask2,2),dataparams.numangles);
k0=zeros(1,dataparams.numangles);

for I=1:dataparams.numangles
    lb=2;
    if dataparams.nrBands==2
        hb=2;
        fb=lb;
    elseif dataparams.nrBands==3
        hb=4;
        fb=hb;
    end

    dataparams.phaOff=0;
    dataparams.fac=ones(1,dataparams.nrBands);
    separateII=separateBands(IIrawFFT(:,:,(I-1)*dataparams.numsteps+1:I*dataparams.numsteps),dataparams.phaOff,dataparams.nrBands,dataparams.fac);

    SeparateII{1,I}=separateII;

    if dataparams.nrBands==2
        c0=separateII(:,:,1);
        c2=separateII(:,:,lb);

        c0=(c0./(max(max(abs(c0)))));
        c2=(c2./(max(max(abs(c2)))));

        c0=c0.*NotchFilter;
        c2=c2.*NotchFilter;

        c0=FFT2D(c0,false);
        c2=FFT2D(c2,false);
        c2=c2.*conj(c0);
        c2=c2./max(max(c2));

        vec=fftshift(FFT2D(c2,true));
    elseif dataparams.nrBands==3
        c0=separateII(:,:,1);
        c3=separateII(:,:,hb);
%         ishow(c0,'0');
%         ishow(separateII(:,:,2),'+1');
%         ishow(separateII(:,:,3),'-1');
%         ishow(c3,'+2');
%         ishow(separateII(:,:,5),'-2');

        c0=c0./(max(max(abs(separateII(:,:,1))))).*NotchFilter;
        c3=c3./(max(max(abs(separateII(:,:,hb))))).*NotchFilter;

        c0=FFT2D(c0,false);
        c3=FFT2D(c3,false);
        c3=c3.*conj(c0);
        c3=c3./max(max(c3));

        vec=fftshift(FFT2D(c3,true));
        ishow(vec,'CrossCorrelation');
    end
    CrossCorrelation(:,:,I)=vec;
    %%
    %         temp=vec.*NotchFilter;
    temp=vec.*NotchFilter2;
    temp=log(1+abs(temp));
    temp=temp./max(max(temp));
    %         MIJ.createImage(temp);

    [yPos,xPos]=find(temp==max(max(temp)));
    peak.xPos=xPos(1);
    peak.yPos=yPos(1);
    k0(I)=sqrt((peak.xPos-cnt(1))^2+(peak.yPos-cnt(2))^2);
    %% 显示findpeak结果
    siz=size(vec);
    w=siz(2);
    h=siz(1);

    x=1:w;
    y=1:h;
    [x,y]=meshgrid(x,y);%x为列,y为行,但矩阵dim=1为行,dim=2为列,k=[peak{1,I}.yPos,peak{1,I}.xPos]
    Rpeak = sqrt( (y-peak.yPos).^2 + (x-peak.xPos).^2 )<5;
    flag=2;
    show_2m=0;
    if(flag==1 || show_2m==1)
        ishow(c1,'c1');   %显示自相关算法后的频谱图片，可看见一个明显亮点
        ishow(Rpeak.*c1,'Fpeak')%显示是否找到的明显亮点
    end

    if(flag==2 || show_2m==1)
        ishow(temp,'temp');   %显示自相关算法后的频谱图片，可看见一个明显亮点
        ishow(Rpeak.*temp,'Fpeak_auto')%显示是否找到的明显亮点
    end
end

Flag=0;
if dataparams.numangles>2           % For very few special cases
    if max(k0)-min(k0)>8
        Flag=1;
        Kobject=min(k0);
        %             Kobject=208;
        Mask1=rad>=(Kobject+1);
        Mask2=rad<=(Kobject-1);
    end
end

for I=1:dataparams.numangles
    vec=CrossCorrelation(:,:,I);
    if Flag==1
        vec(Mask1)=0;
        vec(Mask2)=0;
    end
    temp=vec.*NotchFilter2;
    temp=log(1+abs(temp));
    temp=temp./max(max(temp));
    %         MIJ.createImage(temp);

    [yPos,xPos]=find(temp==max(max(temp)));
    peak.xPos=xPos(1);
    peak.yPos=yPos(1);

    cntrl=zeros(10,30);
    overlap=0.15;
    step=2.5;
    bn1=(dataparams.nrBands-1)*2;
    kx=(peak.xPos-cnt(2));
    ky=(peak.yPos-cnt(1));

    separateII=SeparateII{1,I};
    [peak,cntrl]=fitPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,fb)/(max(max(abs(separateII(:,:,fb))))),1,bn1,OtfProvider,-kx,-ky,overlap,step,cntrl);

    if lb~=hb
        if dataparams.nrBands==2
            peak.kx=peak.kx*2;
            peak.ky=peak.ky*2;
        end

        p1=getPeak(separateII(:,:,1),separateII(:,:,lb),0,1,OtfProvider,peak.kx/2,peak.ky/2,overlap);
        p2=getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,hb)/(max(max(abs(separateII(:,:,1))))),0,2,OtfProvider,peak.kx,peak.ky,overlap);

        dataparams.Dir(I).px=-peak.kx/2;
        dataparams.Dir(I).py=-peak.ky/2;
        dataparams.Dir(I).phaOff=-phase(p1);
        Temp_m1=abs(p1);
        Temp_m2=abs(p2);

        Temp_m1(Temp_m1>1.0)=1;
        Temp_m2(Temp_m2>1.0)=1.0;
        Temp_m1(Temp_m1<0.2)=0.6;
        Temp_m2(Temp_m2<0.2)=0.6;
        dataparams.Dir(I).modul(1)=Temp_m1;
        dataparams.Dir(I).modul(2)=Temp_m2;

        par=dataparams.Dir(I);
        dataparams.fac(2:dataparams.nrBands)=dataparams.Dir(I).modul(1:dataparams.nrBands-1);
        dataparams.fac(2:dataparams.nrBands)=dataparams.Dir(I).modul(1:dataparams.nrBands-1);
        separate_correct=separateBands(IIrawFFT(:,:,(I-1)*dataparams.numsteps+1:I*dataparams.numsteps),par.phaOff,dataparams.nrBands,dataparams.fac);

        SeparateII_cell{1,I} = separate_correct;
    end
    if lb==hb
        p1=getPeak(separateII(:,:,1),separateII(:,:,lb),1,lb,OtfProvider,peak.kx,peak.ky,overlap);
        dataparams.Dir(I).px=-peak.kx;
        dataparams.Dir(I).py=-peak.ky;
        dataparams.Dir(I).phaOff=-phase(p1);
        Temp_m=abs(p1);
        Temp_m(Temp_m>1.0)=1.0;
        Temp_m(Temp_m<0.35)=1-Temp_m;
        dataparams.Dir(I).modul=Temp_m;
        par=dataparams.Dir(I);
        dataparams.fac(2:dataparams.nrBands)=dataparams.Dir(I).modul(1:dataparams.nrBands-1);
        dataparams.fac(2:dataparams.nrBands)=dataparams.Dir(I).modul(1:dataparams.nrBands-1);
        separate_correct=separateBands(IIrawFFT(:,:,(I-1)*dataparams.numsteps+1:I*dataparams.numsteps),par.phaOff,dataparams.nrBands,dataparams.fac);

        SeparateII_cell{1,I} = separate_correct;

    end
    K0(I)=sqrt((dataparams.Dir(I).px)^2+(dataparams.Dir(I).py)^2);
end
%%
for jangle = 1:numangles
    dataparams.allpatternpitch(:,jangle) = [dataparams.Dir(jangle).px,dataparams.Dir(jangle).py];
    dataparams.allpatternangle(:,jangle) = atan(-dataparams.Dir(jangle).py/dataparams.Dir(jangle).px)*180/pi;
    if(numsteps==5)
    dataparams.allpatternphases(:,jangle) = dataparams.Dir(jangle).phaOff + [0 2*pi/5 4*pi/5 6*pi/5 8*pi/5];
    end
    if(numsteps==3)
    dataparams.allpatternphases(:,jangle) = dataparams.Dir(jangle).phaOff + [0 2*pi/3 4*pi/3];
    end
    dataparams.allpatternphases_d(:,jangle) = dataparams.allpatternphases(:,jangle)*180/pi;
    dataparams.allmodule(:,jangle) = [1 dataparams.Dir(jangle).modul(1) dataparams.Dir(jangle).modul(2)];
    
end

SIMparam=zeros(3,5);
if dataparams.numsteps==3
    for i=1:dataparams.numangles
        SIMparam(i,1)=atan(dataparams.Dir(i).py/dataparams.Dir(i).px)*180/pi;
        SIMparam(i,2)=sqrt((dataparams.Dir(i).px)^2+(dataparams.Dir(i).py)^2);
        SIMparam(i,3)=dataparams.Dir(i).phaOff;
        SIMparam(i,4)=dataparams.Dir(i).modul;
    end
elseif dataparams.numsteps==5
    for i=1:dataparams.numangles
        SIMparam(i,1)=atan(dataparams.Dir(i).py/dataparams.Dir(i).px)*180/pi;
        SIMparam(i,2)=sqrt((dataparams.Dir(i).px)^2+(dataparams.Dir(i).py)^2)*2;
        SIMparam(i,3)=dataparams.Dir(i).phaOff;
        SIMparam(i,4)=dataparams.Dir(i).modul(1);
        SIMparam(i,5)=dataparams.Dir(i).modul(2);
    end
end

end