function []=save_as_bin(A,name)
fid=fopen('savedata\'+string(name)+'.bin',"wb");
% writematrix(A,'savedata\'+name+'.txt','Delimiter',' '); 
fwrite(fid,A','float32');
fclose(fid);