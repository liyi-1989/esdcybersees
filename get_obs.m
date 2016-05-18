function [pr,All]=get_obs(mylon,mylat)
    pr=[];
    All=[];
    for i =1979:2003
       fprintf('Year %d\n',i);
       filename=['/media/liyi/新加卷/downscaling/pr_',num2str(i),'.nc'];
       lon=ncread(filename,'lon');
       lat=ncread(filename,'lat'); 
       Pr=ncread(filename,'precipitation_amount');
       tmp=Pr(:,((mylon(1)-0.01)<lon)&(lon<(mylon(2)+0.01)),((mylat(1)-0.01)<lat)&(lat<(mylat(2)+0.01)));
       All=cat(1,All,tmp);
       tmp=mean(tmp,3);
       tmp=mean(tmp,2);
       if length(tmp)==366
           tmp(31+28)=(tmp(31+28)+tmp(31+29))/2;
           tmp(31+29)=[];
       end
       pr=[pr,tmp];
       clear('Pr');
    end
    
    %prm=reshape(pr,365,length(pr)/365);

end