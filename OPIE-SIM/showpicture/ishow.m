function [] = ishow(pic,name)
p = 10;
pic=abs(pic).^(1/p);
minL = min(min(pic));
maxL = max(max( pic ));
figure;
imshow(pic,[minL maxL])
title(name)
