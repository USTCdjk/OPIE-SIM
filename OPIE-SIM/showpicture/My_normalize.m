function OutImg = My_normalize(InImg,pmax,pmin,background_strength)
  if(nargin<2)
      pmax=1;
      pmin=0;
  end
  if(nargin<4)
      background_strength=0.2;
  end
  ymax=pmax;ymin=pmin;
  xmax = max(max(InImg)); %求得InImg中的最大值
  xmin = min(min(InImg)); %求得InImg中的最小值
  OutImg = ((ymax-ymin)*(InImg-xmin)/(xmax-xmin) + ymin)-background_strength*(pmax-pmin); %归一化并取整
  OutImg(OutImg<pmin)=pmin;

  xmax2 = max(max(OutImg)); %求得InImg中的最大值
  xmin2 = min(min(OutImg)); %求得InImg中的最小值
  OutImg = ((ymax-ymin)*(OutImg-xmin2)/(xmax2-xmin2) + ymin); %归一化并取整
end
