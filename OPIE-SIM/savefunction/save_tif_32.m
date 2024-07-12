function [] =save_tif_32(pic,stackfilename)
% maxnum_pic = max(max(pic));
% pic=abs(pic);
final_pic = single(100*pic);
[Nx,Ny,Nz]=size(final_pic);
for j = 1:Nz 
    if j == 1
        t = Tiff(stackfilename,'w');
    else
        t = Tiff(stackfilename, 'a');
    end
    t.setTag('ImageLength', Nx);
    t.setTag('ImageWidth', Ny);
    t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
    t.setTag('BitsPerSample', 32);
    t.setTag('SamplesPerPixel', 1);
    t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
    t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    t.write(final_pic(:,:,j));
    t.close;
end