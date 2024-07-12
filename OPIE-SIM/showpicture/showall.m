function []=showall(data2D)

[~,~,numsteps,Nz,numangles]=size(data2D.allimages_in_Zstack);
for angle=1:numangles
    for j=1:Nz
        for phase=1:numsteps      
            ishow_pic(data2D.allimages_in_Zstack(:,:,phase,j,angle),'phase:'+string(phase)+';Nz:'+string(j)+';angle:'+string(angle)); 
        end
    end
end