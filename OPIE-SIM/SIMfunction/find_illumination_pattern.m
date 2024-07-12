function [dataparams,SeparateII_cell] = find_illumination_pattern(Snoisy,dataparams,ok,OTFflag,OTF_name,attenuation,use_predec,usebackgroundsupression)

%% Basic parameter
numangles = dataparams.numangles;
numsteps = dataparams.numsteps;
numpixelsx = dataparams.numpixelsx;
numpixelsy = dataparams.numpixelsy;
NA = dataparams.NA;
lambda = dataparams.emwavelength;
dataparams.lambda=lambda;
nrBands = dataparams.nrBands;
cyclesPerMicron = double(1./(dataparams.numpixelsx*dataparams.rawpixelsize(1)*0.001));
dataparams.cyclesPerMicron=cyclesPerMicron;
W_opt=0;
%% Calculate OTF
if(ok==0)
    if(OTFflag==1)
        OtfProvider = SimOtfProvider(dataparams,NA,lambda,1);
        PSF_cac=abs(otf2psf((OtfProvider.otf)));
        OTF_cac=OtfProvider.otf;
%         OTF_cac = OTFpost(OTF_cac);
%         OtfProvider.otf=OTF_cac;
        % ishow_3D(OTF_cac,'OTF_cac');
%         ishow(OtfProvider.otf,'Otf')  %显示自动生成的OTF
        % ishow(psf,'psf')
    end
    if(OTFflag==0)
        OtfProvider = SimOtfProvider(dataparams,NA,lambda,cyclesPerMicron,1);
        OTF_cac = double(imread(OTF_name));
        OTF_cac = OTFpost(OTF_cac);
        OtfProvider.otf=OTF_cac;
        PSF_cac=abs(otf2psf((OTF_cac)));
        ishow_3D(OTF_cac,'OTF_cac');
        % ishow(OtfProvider.otf,'Otf')  %显示OTF
        % ishow(psf,'psf')
    end
    dataparams.OtfProvider=OtfProvider;
    dataparams.OTF_cac=OTF_cac;
%     save_as_mat(OTF_cac,'OTF_cac');
    dataparams.PSF_cac=PSF_cac;
end
if(ok~=0)
    OtfProvider=dataparams.OtfProvider;
    PSF_cac=dataparams.PSF_cac;
end

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
cnt=[numpixelsx/2+1,numpixelsx/2+1];
cutoff=1000/(0.5*lambda/NA);
[x,y]=meshgrid(1:numpixelsx,1:numpixelsy);
rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);
Mask1=double(rad<=1.0*(cutoff/cyclesPerMicron+1));
NotchFilter0=getotfAtt(numpixelsx,cyclesPerMicron,0.3*cutoff,0,0);
% ishow(NotchFilter0,'NF0');%可视化滤波器

NotchFilter1=NotchFilter0.*Mask1;
% ishow(NotchFilter1,'NF1');

Mask2=double(rad<=1.10*(cutoff/cyclesPerMicron+1));

NotchFilter2=NotchFilter0.*Mask2;
% ishow(NotchFilter2,'NF2');

%% Cross-correlation to solve illumination pattern
% CrossCorrelation=zeros(size(Mask2,1),size(Mask2,2),numangles);
% k0=zeros(1,numangles); 

for I=1:numangles
    cntrl = zeros(10,30);
    overlap = 0.15;
    step = 2.5;
    bn1 = (nrBands-1)*2;
%     ishow(IIrawFFT(:,:,(I-1)*numsteps+1),'o');
%% 频谱分离，见北大第二页III. SIM FORMULATION
    separateII = SeparatedComponents2D(IIrawFFT(:,:,(I-1)*numsteps+1:I*numsteps),numsteps,0,0,1);
    SeparateII_cell{1,I} = separateII; %separateII为北大III. SIM FORMULATION公式（4）里的S(k) H(k)，S(k − pθ) H(k)，S(k + pθ) H(k)
    
    %% 相关算法找peak
    c0 = separateII(:,:,1);
    c1 = separateII(:,:,2);
    c2 = separateII(:,:,3);

%     save_as_mat(c0,'seperate_fo_dir'+string(I));
%     save_as_mat(c1,'seperate_fp_dir'+string(I));
%     save_as_mat(c2,'seperate_fm_dir'+string(I));

    ishow(c0,'o'); %显示频谱分离后结果，o应该只有中心亮点，p与m在两侧有两点
    ishow(c1,'p');
    ishow(c2,'m');

    c0 = c0./(max(max(abs(separateII(:,:,1))))).*NotchFilter2;
    c1 = c1./(max(max(abs(separateII(:,:,2))))).*NotchFilter2;

%%直接找peak
    [yPos_1,xPos_1] = find(c1==max(max(c1)));
    peak{1,I}.xPos_1 = xPos_1(1);
    peak{1,I}.yPos_1 = yPos_1(1);

    %     c2 = c2./(max(max(abs(separateII(:,:,3))))).*NotchFilter1;
%     ishow(c0,'0q');
%     ishow(c1,'1q');
    %     ishow(c2,'2q');

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
    temp = temp./max(max(temp));

    ishow(temp,'temp');
    [yPos_2,xPos_2] = find(temp==max(max(temp)));
    peak{1,I}.xPos_2 = xPos_2(1);
    peak{1,I}.yPos_2 = yPos_2(1);
%     k0(I) = sqrt((peak{1,I}.xPos-cnt(1))^2+(peak{1,I}.yPos-cnt(2))^2);
end


for I = 1:numangles
    a=peak{1,I};
    kx_1 = (a.xPos_1-cnt(2));
    ky_1 = (a.yPos_1-cnt(1));
    freq_1(:,I) = [-kx_1,-ky_1];
    ang_1(:,I) = atan(ky_1/kx_1)*180/pi;
end

evaluate_freq1_1 = [sqrt(freq_1(1,1)^2+freq_1(2,1)^2),sqrt(freq_1(1,2)^2+freq_1(2,2)^2),sqrt(freq_1(1,3)^2+freq_1(2,3)^2)];
evaluate_freq_std1_1 = std(evaluate_freq1_1,1);
evaluate_ang1_1 = [mod(ang_1(1,1)-ang_1(1,2),180),mod(ang_1(1,2)-ang_1(1,3),180),mod(ang_1(1,3)-ang_1(1,1),180)];
evaluate_ang_std1_1 = std(evaluate_ang1_1,1);

for I = 1:numangles
    kx_2 = (peak{1,I}.xPos_2-cnt(2));
    ky_2 = (peak{1,I}.yPos_2-cnt(1));
    freq_2(:,I) = [-kx_2,-ky_2];
    ang_2(:,I) = atan(ky_2/kx_2)*180/pi;
end

evaluate_freq1_2 = [sqrt(freq_2(1,1)^2+freq_2(2,1)^2),sqrt(freq_2(1,2)^2+freq_2(2,2)^2),sqrt(freq_2(1,3)^2+freq_2(2,3)^2)];
evaluate_freq_std1_2 = std(evaluate_freq1_2,1);
evaluate_ang1_2 = [mod(ang_2(1,1)-ang_2(1,2),180),mod(ang_2(1,2)-ang_2(1,3),180),mod(ang_2(1,3)-ang_2(1,1),180)];
evaluate_ang_std1_2 = std(evaluate_ang1_2,1);

if (evaluate_freq_std1_1>0.5 && evaluate_ang_std1_1>5)&&(evaluate_freq_std1_2>0.5 && evaluate_ang_std1_2>5)
    disp("fail to find frequency, phase and modulation depth");
end
if evaluate_freq_std1_2<0.5||evaluate_ang_std1_2<5
    for I = 1:numangles
        peak{1,I}.xPos = peak{1,I}.xPos_2(1);
        peak{1,I}.yPos = peak{1,I}.yPos_2(1);
        flag=2;
    end
elseif evaluate_freq_std1_1<0.5||evaluate_ang_std1_1<5
    for I = 1:numangles
        peak{1,I}.xPos = peak{1,I}.xPos_1(1);
        peak{1,I}.yPos = peak{1,I}.yPos_1(1);
        flag=1;
    end
end

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
    %% 显示findpeak结果
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
    [peak{1,I},~] = fitPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,2)/(max(max(abs(separateII(:,:,2))))),1,bn1,OtfProvider,-kx,-ky,overlap,step,cntrl);
end
tic
for I=1:numangles
%     ishow(IIrawFFT(:,:,(I-1)*numsteps+1),'o');
%% 频谱分离，见北大第二页III. SIM FORMULATION
    separateII=SeparateII_cell{1,I} ;    

    if(ok~=0)
        peak{1,I}=dataparams.peak_initial{1,I};
    end
%     ishow(FFT2D(separateII(:,:,1),true),'separateII(:,:,1)');
%%  相关算法计算phaseoff和调制度，参考北大第五页G.phase matching
    p1 = getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,2)/(max(max(abs(separateII(:,:,1))))),0,1,OtfProvider,peak{1,I}.kx,peak{1,I}.ky,overlap);
%     p1 = getPeak(separateII(:,:,1),separateII(:,:,2),0,1,OtfProvider,peak{1,I}.kx,peak{1,I}.ky,overlap,I);
    %前两个参数对应公式（16）Su(k)与Ss(k − pθ)
    
    params.Dir(I).px = -peak{1,I}.kx;
    params.Dir(I).py = -peak{1,I}.ky;

%% 迭代计算phase
    params.k_OPT(:,I)=[peak{1,I}.ky,peak{1,I}.kx]';
    for nsteps=1:numsteps
        [phase_OPT] = IlluminationPhaseF(IIraw(:,:,nsteps+numsteps*(I-1)),params.k_OPT(:,I)');
        params.phase_OPT(nsteps,I)=phase_OPT;
    end
    
    % for display in command window
    phase_OPT_d= params.phase_OPT(:,I)*180/pi;
    params.phase_OPT_d(:,I)=phase_OPT_d;

    if(numsteps==5)
    evaluate_phaseoff=mod((params.phase_OPT(:,I)'-[params.phase_OPT(5,I) params.phase_OPT(1,I) params.phase_OPT(2,I) params.phase_OPT(3,I) params.phase_OPT(4,I)]),2*pi);
    end
    if(numsteps==3)
    evaluate_phaseoff=mod((params.phase_OPT(:,I)'-[0 2*pi/3 4*pi/3]),2*pi);
    end

    evaluate_phaseoff_d=evaluate_phaseoff*180/pi;
    disp(evaluate_phaseoff_d);
    evaluate_phaseoff_std1 = std(evaluate_phaseoff,1);
    if(evaluate_phaseoff_std1>5 && ok==0)
        disp('phaseoff caculate error')
    end
    
    if(W_opt==1)
        params.Dir(I).phaOff=mean(evaluate_phaseoff);
    end
%%  相关算法计算phaseoff和调制度，参考北大第五页G.phase matching
    if(W_opt==0)
        params.Dir(I).phaOff=-angle(p1);
        disp(-angle(p1));
    end
    Temp_m1 = abs(p1);
    
    if(Temp_m1>1.0)
        Temp_m1= 1.0;
%     elseif(Temp_m1<0.5 && Temp_m1>0.2)
%         Temp_m1 =0.5;
    elseif(Temp_m1<0.15)
        Temp_m1 =0.6;
    end
    params.Dir(I).modul(1) = Temp_m1;
    
    separateII(:,:,2)=exp(+1i*params.Dir(I).phaOff)./Temp_m1.*separateII(:,:,2);
    separateII(:,:,3)=exp(-1i*params.Dir(I).phaOff)./Temp_m1.*separateII(:,:,3);
    SeparateII_cell{1,I} = separateII;

    
%     separateII = SeparateII_cell{1,I};
%     p2 = getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,2)/(max(max(abs(separateII(:,:,1))))),0,1,OtfProvider,peak{1,I}.kx,peak{1,I}.ky,overlap);
%     p3 = getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,3)/(max(max(abs(separateII(:,:,1))))),0,1,OtfProvider,-peak{1,I}.kx,-peak{1,I}.ky,overlap);
%     disp(-angle(p2));%消除phaseoff后理论上p2，p3是实数，显示的值应该非常小
%     disp(-angle(p3));
    disp('频谱分离:'+string(I))
end
toc

%% Save the results
for jangle = 1:numangles
    dataparams.allpatternpitch(:,jangle) = [params.Dir(jangle).px,params.Dir(jangle).py];
    dataparams.Dir=params.Dir;
    
    dataparams.allpatternangle(:,jangle) = atan(-params.Dir(jangle).py/params.Dir(jangle).px)*180/pi;
    if(numsteps==5)
    dataparams.allpatternphases(:,jangle) = params.Dir(jangle).phaOff + [0 2*pi/5 4*pi/5 6*pi/5 8*pi/5];
    end
    if(numsteps==3)
    dataparams.allpatternphases(:,jangle) = params.Dir(jangle).phaOff + [0 2*pi/3 4*pi/3];
    end
    dataparams.allpatternphases_d(:,jangle) = dataparams.allpatternphases(:,jangle)*180/pi;
    dataparams.allmodule(:,jangle) = [1 params.Dir(jangle).modul(1)];
end

for I = 1:numangles
    if(W_opt==1)
        separateII = SeparatedComponents2D(IIrawFFT(:,:,(I-1)*numsteps+1:I*numsteps),numsteps,1,params.phase_OPT(:,I),params.Dir(I).modul(1));
        SeparateII_cell{1,I} = separateII; %separateII为北大III. SIM FORMULATION公式（4）里的S(k) H(k)，S(k − pθ) H(k)，S(k + pθ) H(k)
    end
%     if(W_opt==0)
%         separateII = SeparatedComponents2D(IIrawFFT(:,:,(I-1)*numsteps+1:I*numsteps),numsteps,1,dataparams.allpatternphases(:,I),params.Dir(I).modul(1));
%         SeparateII_cell{1,I} = separateII; %separateII为北大III. SIM FORMULATION公式（4）里的S(k) H(k)，S(k − pθ) H(k)，S(k + pθ) H(k)
%     end
end
dataparams.lambda = dataparams.exwavelength;
dataparams.OtfProvider = SimOtfProvider(dataparams,NA,dataparams.lambda,attenuation);
% dataparams.OtfProvider=OtfProvider;
end