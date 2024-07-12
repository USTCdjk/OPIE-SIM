function []=save_as_mat(A,name)

save('savedata\'+string(name)+'.mat','A')