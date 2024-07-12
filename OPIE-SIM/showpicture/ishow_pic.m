function [] = ishow_pic(pic,name)
p = 1;
pic=abs(pic).^(1/p);
minL = min(min(pic));
maxL = max(max( pic ));
figure;
imshow(pic,[minL maxL])
title(name)