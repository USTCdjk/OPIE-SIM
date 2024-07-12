function [] =save_tif_app(pic,stackfilename)
maxnum_pic = max(max(pic));
final_pic = uint8(255*pic./maxnum_pic);
imwrite(final_pic, stackfilename,'WriteMode','append') % 写入stack图像