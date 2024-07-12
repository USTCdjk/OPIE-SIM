function [final_pic] = generateresult(data2D,usebackgroundsupression)
final_pic=sum(data2D.final_image(:,:,:),3);
% final_pic(final_pic<0.01*min(min(final_pic)))=0;
if(usebackgroundsupression)
    final_pic=My_normalize(final_pic);
end
save_tif_32(final_pic,'.\save\OTFcorrection_final_image.tif')
end