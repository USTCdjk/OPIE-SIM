function FWHM = PSF_FWHM(PSFo)

w = size(PSFo,1);
FWHM_1=0;
FWHM_2=0;

PSFo(PSFo<0.5*max(max(PSFo)))=0;
PSF_1 = sum(PSFo,1);
PSF_2 = sum(PSFo,2);
for i=1:w
    if ( abs(PSF_1(1,i))>0 )
    	FWHM_1 = FWHM_1+1;
    end
    if ( abs(PSF_2(i,1))>0 )
    	FWHM_2 = FWHM_2+1;
    end
end
FWHM=max(FWHM_1,FWHM_2);