function [separate] = SeparatedComponents2D( IrawFFT,numsteps,corr,corr_phases,MF)
% Aim: Unmixing the frequency components of raw SIM images
%   phaseShift,phaseShift0: illumination phase shifts
%   FcS1aT,FcS2aT,FcS3aT: FT of raw SIM images
%   FiSMao,FiSMap,FiSMam: unmixed frequency components of raw SIM images
if(corr==0)
phaPerBand=numsteps;
phases=zeros(1,phaPerBand);
for p=1:phaPerBand
    phases(p)=(2*pi*(p-1))/phaPerBand;%用来生成北大III. SIM FORMULATION公式（4）里的M，此处phaOff=0，没用准确的相位，具体影响参考"结构光算法/《相位影响推导》"
end
end
if(corr==1)
    phases=corr_phases;
end
if(numsteps==3)
    phaseShift0 = phases(1);
    phaseShift1 = phases(2);
    phaseShift2 = phases(3);

    phase_d=[phaseShift0*180/pi,phaseShift1*180/pi,phaseShift2*180/pi];
    %% Transformation Matrix
    M = 0.5*[1 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0);
        1 0.5*MF*exp(-1i*phaseShift1) 0.5*MF*exp(+1i*phaseShift1);
        1 0.5*MF*exp(-1i*phaseShift2) 0.5*MF*exp(+1i*phaseShift2)];%北大III. SIM FORMULATION公式（4）里的M,3相位情况

    %% Separting the components
    %===========================================================
    Minv = inv(M);%求逆矩阵

    FiSMao = Minv(1,1)*IrawFFT(:,:,1) + Minv(1,2)*IrawFFT(:,:,2) + Minv(1,3)*IrawFFT(:,:,3);
    FiSMap = Minv(2,1)*IrawFFT(:,:,1) + Minv(2,2)*IrawFFT(:,:,2) + Minv(2,3)*IrawFFT(:,:,3);
    FiSMam = Minv(3,1)*IrawFFT(:,:,1) + Minv(3,2)*IrawFFT(:,:,2) + Minv(3,3)*IrawFFT(:,:,3);
    
    separate(:,:,1)=FiSMao;
    separate(:,:,2)=FiSMap;
    separate(:,:,3)=FiSMam;
end

if(numsteps==5)
    phaseShift0 = phases(1);
    phaseShift1 = phases(2);
    phaseShift2 = phases(3);
    phaseShift3 = phases(4);
    phaseShift4 = phases(5);

    phase_d=[phaseShift0*180/pi,phaseShift1*180/pi,phaseShift2*180/pi,phaseShift3*180/pi,phaseShift4*180/pi];
    %% Transformation Matrix
    M = 0.5*[1 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0);
        1 0.5*MF*exp(-1i*phaseShift1) 0.5*MF*exp(+1i*phaseShift1);
        1 0.5*MF*exp(-1i*phaseShift2) 0.5*MF*exp(+1i*phaseShift2);
        1 0.5*MF*exp(-1i*phaseShift3) 0.5*MF*exp(+1i*phaseShift3);
        1 0.5*MF*exp(-1i*phaseShift4) 0.5*MF*exp(+1i*phaseShift4)];%北大III. SIM FORMULATION公式（4）里的M,5相位情况

    %% Separting the components
    %===========================================================
    M_temp1 = M'*M;
    Minv = M_temp1\M';%伪逆矩阵，利用计算方法中的极大似然法
    
%     save_as_mat(Minv,'M^-1');

    FiSMao = Minv(1,1)*IrawFFT(:,:,1) + Minv(1,2)*IrawFFT(:,:,2) + Minv(1,3)*IrawFFT(:,:,3) + Minv(1,4)*IrawFFT(:,:,4) + Minv(1,5)*IrawFFT(:,:,5);
    FiSMap = Minv(2,1)*IrawFFT(:,:,1) + Minv(2,2)*IrawFFT(:,:,2) + Minv(2,3)*IrawFFT(:,:,3) + Minv(2,4)*IrawFFT(:,:,4) + Minv(2,5)*IrawFFT(:,:,5);
    FiSMam = Minv(3,1)*IrawFFT(:,:,1) + Minv(3,2)*IrawFFT(:,:,2) + Minv(3,3)*IrawFFT(:,:,3) + Minv(3,4)*IrawFFT(:,:,4) + Minv(3,5)*IrawFFT(:,:,5);
    
    separate(:,:,1)=FiSMao;
    separate(:,:,2)=FiSMap;
    separate(:,:,3)=FiSMam;
end

