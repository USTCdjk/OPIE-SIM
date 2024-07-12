function [] = ishow_Re(pic,name)
pic = real( ifft2(fftshift(pic)) );
pic =My_normalize(pic,3,0);
minL = min(min(pic));
maxL = max(max( pic ));
figure;
imshow(pic,[minL maxL])
title(name)