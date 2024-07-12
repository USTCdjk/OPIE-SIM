function []=save_as_bin_1b(A,name)
filename='.\savedata\'+string(name)+'.bin';
fid=fopen(filename,'wb');
if fid==-1
    disp('打开失败');
end
% writematrix(A,'savedata\'+string(name)+'.txt','Delimiter',' '); 
fwrite(fid,A','ubit1');
fclose(fid);