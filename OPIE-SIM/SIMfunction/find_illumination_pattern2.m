function [dataparams,SeparateII_cell,find_peak_ok] = find_illumination_pattern2(Snoisy,dataparams,ok,OTFflag,OTF_name,attenuation,use_predec,usebackgroundsupression)

%% Basic parameter
numangles = dataparams.numangles;
numsteps = dataparams.numsteps;
numpixelsx = dataparams.numpixelsx;
numpixelsy = dataparams.numpixelsy;
NA = dataparams.NA;
lambda = dataparams.emwavelength;
dataparams.lambda = lambda;
nrBands = dataparams.nrBands;
cyclesPerMicron = double(1./(dataparams.numpixelsx*dataparams.rawpixelsize(1)*0.001));
dataparams.cyclesPerMicron =cyclesPerMicron;
W_opt=1;
%% Calculate OTF

OtfProvider = SimOtfProvider(dataparams,NA,lambda,1);
save_as_mat(OtfProvider.otf,'OTFideal')
PSF_cac=abs(otf2psf((OtfProvider.otf)));
save_tif_32(OtfProvider.otf,'.\save\Otf.tif');
save_tif_32(PSF_cac,'.\save\PSF_cac.tif');

dataparams.OtfProvider=OtfProvider;

if(use_predec==1)
    Temp=importImages(Snoisy);
    IIraw=deconvlucy(Temp,PSF_cac,5);
else
    IIraw=Snoisy;
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

%% Cutoff and notchfilter
cnt=[numpixelsx/2+1,numpixelsy/2+1];
cutoff=1000/(0.5*lambda/NA);
[x,y]=meshgrid(1:numpixelsx,1:numpixelsy);
radx=(x-cnt(2)).^2;
rady=(y-cnt(1)).^2;
maskxy=double(radx>-1).*double(rady>-1);

rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);
Mask1=double(rad<=1.0*(cutoff/cyclesPerMicron+1));
NotchFilter01=getotfAtt(numpixelsx,cyclesPerMicron,1.0*cutoff,0,0);
NotchFilter02=getotfAtt(numpixelsx,cyclesPerMicron,0.5*cutoff,0,0);
% NotchFilter03=getotfAtt(numpixelsx,10*cyclesPerMicron,0.05*cutoff,0,0);
NotchFilter1=NotchFilter01.*Mask1.*maskxy;
Mask2=double(rad<=1.10*(cutoff/cyclesPerMicron+1));
NotchFilter2=NotchFilter02.*Mask2.*maskxy;

%% Cross-correlation to solve illumination pattern
% CrossCorrelation=zeros(size(Mask2,1),size(Mask2,2),numangles);
% k0=zeros(1,numangles);
cntrl = zeros(10,30);
overlap = 0.15;
step = 2.5;
bn1 = (nrBands-1)*2;

%% 迭代计算phase
if(W_opt==1)
    for I=1:numangles
        if(numsteps==5)
        c1=IIrawFFT(:,:,(I-1)*numsteps+1);
%         c_show=c1.*NotchFilter03;
        c1 = c1./(max(max(abs(c1)))).*NotchFilter1;
        c2=IIrawFFT(:,:,(I-1)*numsteps+2);
%         c_show=c_show.*c2.*NotchFilter03;
        c2 = c2./(max(max(abs(c2)))).*NotchFilter1;
        c3=IIrawFFT(:,:,(I-1)*numsteps+3);
%         c_show=c_show.*c3.*NotchFilter03;
        c3 = c3./(max(max(abs(c3)))).*NotchFilter1;
        c4=IIrawFFT(:,:,(I-1)*numsteps+4);
%         c_show=c_show.*c4.*NotchFilter03;
        c4 = c4./(max(max(abs(c4)))).*NotchFilter1;
        c5=IIrawFFT(:,:,(I-1)*numsteps+5);
%         c_show=c_show.*c5.*NotchFilter03;
        c5 = c5./(max(max(abs(c5)))).*NotchFilter1;
        c=c1.*c2.*c3.*c4.*c5;
        
        end
        if(numsteps==3)
        c1=IIrawFFT(:,:,(I-1)*numsteps+1);
%         c_show=c1.*NotchFilter03;
        c1 = c1./(max(max(abs(c1)))).*NotchFilter1;
        c2=IIrawFFT(:,:,(I-1)*numsteps+2);
%         c_show=c_show.*c2.*NotchFilter03;
        c2 = c2./(max(max(abs(c2)))).*NotchFilter1;
        c3=IIrawFFT(:,:,(I-1)*numsteps+3);
%         c_show=c_show.*c3.*NotchFilter03;
        c3 = c3./(max(max(abs(c3)))).*NotchFilter1;
        c=c1.*c2.*c3;
        end
        c = c./(max(max(abs(c))));
        %%直接找peak
        [yPos,xPos] = find(c==max(max(c)));
        peak_opt{1,I}.xPos = xPos(1);
        peak_opt{1,I}.yPos = yPos(1);
        
        siz=size(c);
        w=siz(2);
        h=siz(1);

        x=1:w;
        y=1:h;
        [x,y]=meshgrid(x,y);%x为列,y为行,但矩阵dim=1为行,dim=2为列,k=[peak{1,I}.yPos,peak{1,I}.xPos]
        Rpeak1 = double(sqrt( (y-peak_opt{1,I}.yPos).^2 + (x-peak_opt{1,I}.xPos).^2 )<5);
        Rpeak2 = double(sqrt( (y-h-1+peak_opt{1,I}.yPos).^2 + (x-w-1+peak_opt{1,I}.xPos).^2 )<5);
        ishow(c,'c');   %显示自相关算法后的频谱图片，可看见一个明显亮点
        ishow(Rpeak1.*c+Rpeak2.*c,'Fpeak')%显示是否找到的明显亮点
    end
    k_OPT=zeros(2,numangles);
    for I=1:numangles
        kx = -(peak_opt{1,I}.xPos-cnt(2));
        ky = -(peak_opt{1,I}.yPos-cnt(1));
        k_OPT(:,I)=[ky,kx]';
        for nsteps=1:numsteps
            [phase_OPT] = IlluminationPhaseF(IIraw(:,:,nsteps+numsteps*(I-1)),k_OPT(:,I)');
            params.phase_OPT(nsteps,I)=phase_OPT;
            % for display in command window
            params.phase_OPT_d(nsteps,I)= phase_OPT*180/pi;
        end    
    end
    for I = 1:numangles
        separateII = SeparatedComponents2D(IIrawFFT(:,:,(I-1)*numsteps+1:I*numsteps),numsteps,1,params.phase_OPT(:,I),1);
        SeparateII_cell{1,I} = separateII; %separateII为北大III. SIM FORMULATION公式（4）里的S(k) H(k)，S(k − pθ) H(k)，S(k + pθ) H(k)
    end
    for I = 1:numangles
        dataparams.allpatternpitch(:,I) = [-k_OPT(2,I),-k_OPT(1,I)];
        dataparams.allpatternangle(:,I) = atan(-k_OPT(1,I)/k_OPT(2,I))*180/pi;
        dataparams.allpatternphases(:,I) = params.phase_OPT(:,I);
        dataparams.allmodule(:,I) = [1 0.5];
        if(numsteps==5)
            evaluate_phaseoff=mod((params.phase_OPT(:,I)'-[params.phase_OPT(5,I) params.phase_OPT(1,I) params.phase_OPT(2,I) params.phase_OPT(3,I) params.phase_OPT(4,I)]),2*pi);
        end
        if(numsteps==3)
            evaluate_phaseoff=mod((params.phase_OPT(:,I)'-[params.phase_OPT(3,I) params.phase_OPT(1,I) params.phase_OPT(2,I)]),2*pi);
        end
        evaluate_phaseoff_d=evaluate_phaseoff*180/pi;
        for i=1:numsteps
            if(evaluate_phaseoff_d(i)>180)
                evaluate_phaseoff_d(i)=360-evaluate_phaseoff_d(i);
            end
        end
        disp(evaluate_phaseoff_d);
    end
end

if(W_opt==0)
    for I=1:numangles
        %% 频谱分离，见北大第二页III. SIM FORMULATION
        separateII = SeparatedComponents2D(IIrawFFT(:,:,(I-1)*numsteps+1:I*numsteps),numsteps,0,0,1);
        SeparateII_cell{1,I} = separateII; %separateII为北大III. SIM FORMULATION公式（4）里的S(k) H(k)，S(k − pθ) H(k)，S(k + pθ) H(k)
    end
end

for I=1:numangles
    separateII=SeparateII_cell{1,I} ;
    %% 相关算法找peak
    c0 = separateII(:,:,1);
    c1 = separateII(:,:,2);
    %     c2 = separateII(:,:,3);

    c0 = c0./(max(max(abs(separateII(:,:,1))))).*NotchFilter2;
    c1 = c1./(max(max(abs(separateII(:,:,2))))).*NotchFilter2;


    c0_i = FFT2D(c0,false);
    c1_i = FFT2D(c1,false);
    %     c2_i = FFT2D(c2,false);
    %     ishow(c0,'0');
    %     ishow(c3,'3');

    c0_1 = c1_i.*conj(c0_i);
    c0_1 = c0_1./max(max(c0_1));
    %     ishow(c0_1,'0_1');

    vec = fftshift(FFT2D(c0_1,true));
    %     ishow(vec,'v');

    %     CrossCorrelation(:,:,I) = vec;
    temp = vec.*NotchFilter1;
    temp = log(1+abs(temp));
    %     temp = temp./max(max(temp));


    [yPos_2,xPos_2] = find(temp==max(max(temp)));
    peak{1,I}.xPos_2 = xPos_2(1);
    peak{1,I}.yPos_2 = yPos_2(1);
end

%% 找合适的波矢量
for I = 1:numangles
    freq_1(:,I) = k_OPT(:,I);
    ang_1(:,I) = atan(k_OPT(1,I)/k_OPT(2,I))*180/pi;
end

evaluate_freq1_1 = [sqrt(freq_1(1,1)^2+freq_1(2,1)^2),sqrt(freq_1(1,2)^2+freq_1(2,2)^2),sqrt(freq_1(1,3)^2+freq_1(2,3)^2)];
evaluate_freq_std1_1 = std(evaluate_freq1_1,1);
evaluate_ang1_1 = [mod(ang_1(1,1)-ang_1(1,2),180),mod(ang_1(1,2)-ang_1(1,3),180),mod(ang_1(1,3)-ang_1(1,1),180)];
for i=1:numangles
    if(evaluate_ang1_1(i)>100)
        evaluate_ang1_1(i)=evaluate_ang1_1(i)-60;
    end
end
evaluate_ang_std1_1 = std(evaluate_ang1_1,1);

for I = 1:numangles
    kx_2 = (peak{1,I}.xPos_2-cnt(2));
    ky_2 = (peak{1,I}.yPos_2-cnt(1));
    freq_2(:,I) = [-ky_2,-kx_2];
    ang_2(:,I) = atan(ky_2/kx_2)*180/pi;
end

evaluate_freq1_2 = [sqrt(freq_2(1,1)^2+freq_2(2,1)^2),sqrt(freq_2(1,2)^2+freq_2(2,2)^2),sqrt(freq_2(1,3)^2+freq_2(2,3)^2)];
evaluate_freq_std1_2 = std(evaluate_freq1_2,1);
evaluate_ang1_2 = [mod(ang_2(1,1)-ang_2(1,2),180),mod(ang_2(1,2)-ang_2(1,3),180),mod(ang_2(1,3)-ang_2(1,1),180)];
for i=1:numangles
    if(evaluate_ang1_2(i)>100)
        evaluate_ang1_2(i)=evaluate_ang1_2(i)-60;
    end
end
evaluate_ang_std1_2 = std(evaluate_ang1_2,1);

if (evaluate_freq_std1_1>0.5 && evaluate_ang_std1_1>5)&&(evaluate_freq_std1_2>0.5 && evaluate_ang_std1_2>5)
    disp("fail to find frequency, phase and modulation depth");
    find_peak_ok=0;
end
if evaluate_freq_std1_2<0.5||evaluate_ang_std1_2<5
    for I = 1:numangles
        peak{1,I}.xPos = peak{1,I}.xPos_2(1);
        peak{1,I}.yPos = peak{1,I}.yPos_2(1);
        flag=2;
        find_peak_ok=1;
    end
end
if evaluate_freq_std1_1<0.5||evaluate_ang_std1_1<5||evaluate_freq_std1_2>0.5||evaluate_ang_std1_2>5
    for I = 1:numangles
        peak{1,I}.xPos = peak_opt{1,I}.xPos(1);
        peak{1,I}.yPos = peak_opt{1,I}.yPos(1);
        flag=1;
        find_peak_ok=1;
    end
end
%% 显示findpeak结果
for I=1:numangles
    separateII=SeparateII_cell{1,I} ;
    c0 = separateII(:,:,1);
    c1 = separateII(:,:,2);


    c0 = c0./(max(max(abs(separateII(:,:,1))))).*NotchFilter2;
    c1 = c1./(max(max(abs(separateII(:,:,2))))).*NotchFilter2;


    c0_i = FFT2D(c0,false);
    c1_i = FFT2D(c1,false);


    c0_1 = c1_i.*conj(c0_i);
    c0_1 = c0_1./max(max(c0_1));

    vec = fftshift(FFT2D(c0_1,true));

    %     CrossCorrelation(:,:,I) = vec;
    temp = vec.*NotchFilter1;
    temp = log(1+abs(temp));
    %     temp = temp./max(max(temp));

    siz=size(vec);
    w=siz(2);
    h=siz(1);

    x=1:w;
    y=1:h;
    [x,y]=meshgrid(x,y);%x为列,y为行,但矩阵dim=1为行,dim=2为列,k=[peak{1,I}.yPos,peak{1,I}.xPos]
    Rpeak = sqrt( (y-peak{1,I}.yPos).^2 + (x-peak{1,I}.xPos).^2 )<5;

    show_2m=0;
    if(flag==1 || show_2m==1)
        ishow(c1,'c1');   %显示自相关算法后的频谱图片，可看见一个明显亮点
        ishow(Rpeak.*c1,'Fpeak')%显示是否找到的明显亮点
    end

    if(flag==2 || show_2m==1)
        ishow(temp,'temp');   %显示自相关算法后的频谱图片，可看见一个明显亮点
        ishow(Rpeak.*temp,'Fpeak_auto')%显示是否找到的明显亮点
    end
    %% 计算phaseoff和调制度m

    kx = (peak{1,I}.xPos-cnt(2));
    ky = (peak{1,I}.yPos-cnt(1));

    separateII = SeparateII_cell{1,I};

    %%  迭代3次找亚像素级别的波矢量值
    [peak{1,I},~] = fitPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,2)/(max(max(abs(separateII(:,:,2))))),1,bn1,dataparams.OtfProvider,-kx,-ky,overlap,step,cntrl);
end

for I=1:numangles
    %% 频谱分离，见北大第二页III. SIM FORMULATION
    separateII=SeparateII_cell{1,I} ;
    %%  相关算法计算phaseoff和调制度，参考北大第五页G.phase matching
    %     p1 = getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,2)/(max(max(abs(separateII(:,:,1))))),0,1,OtfProvider,peak{1,I}.kx,peak{1,I}.ky,overlap);
    p1 = getPeak(separateII(:,:,1),separateII(:,:,2),0,1,dataparams.OtfProvider,peak{1,I}.kx,peak{1,I}.ky,overlap);
    %前两个参数对应公式（16）Su(k)与Ss(k − pθ)

    params.Dir(I).px = -peak{1,I}.kx;
    params.Dir(I).py = -peak{1,I}.ky;



    %%  相关算法计算phaseoff和调制度，参考北大第五页G.phase matching

    params.Dir(I).phaOff=-angle(p1);

    Temp_m1 = abs(p1);

    Temp_m1(Temp_m1>1.0) = 1.0;
    Temp_m1(Temp_m1<0) = 0;
    Temp_m1=1-Temp_m1;
    params.Dir(I).modul(1) = Temp_m1;

    separateII(:,:,2)=(exp(+1i*params.Dir(I).phaOff)/Temp_m1).*separateII(:,:,2);
    separateII(:,:,3)=(exp(-1i*params.Dir(I).phaOff)/Temp_m1).*separateII(:,:,3);
    SeparateII_cell{1,I} = separateII;
    ishow(separateII(:,:,1),'o');
    ishow(separateII(:,:,2),'p');
    ishow(separateII(:,:,3),'m');

    %     separateII = SeparateII_cell{1,I};
    %     p2 = getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,2)/(max(max(abs(separateII(:,:,1))))),0,1,OtfProvider,peak{1,I}.kx,peak{1,I}.ky,overlap);
    %     p3 = getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,3)/(max(max(abs(separateII(:,:,1))))),0,1,OtfProvider,-peak{1,I}.kx,-peak{1,I}.ky,overlap);
    %     disp(-angle(p2));%消除phaseoff后理论上p2，p3是实数，显示的值应该非常小
    %     disp(-angle(p3));
    disp('频谱分离:'+string(I))
end


%% Save the results
for jangle = 1:numangles
    dataparams.Dir=params.Dir;
    dataparams.allpatternpitch(:,jangle) = [params.Dir(jangle).px,params.Dir(jangle).py];
    dataparams.allpatternangle(:,jangle) = atan(-params.Dir(jangle).py/params.Dir(jangle).px)*180/pi;
    dataparams.allpatternphases(:,jangle) = params.Dir(jangle).phaOff + dataparams.allpatternphases(:,jangle);
    dataparams.allmodule(:,jangle) = [1 params.Dir(jangle).modul(1)];
end

% dataparams.lambda = dataparams.exwavelength;
% dataparams.OtfProvider = SimOtfProvider(dataparams,NA,dataparams.lambda,attenuation);
% dataparams.OtfProvider=OtfProvider;
end
