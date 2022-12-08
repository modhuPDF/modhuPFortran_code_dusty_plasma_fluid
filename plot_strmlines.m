
%     contour(x,z,cont_shi_d,5)   


% return
 cmax = max(cont_shi_d(:));
 cmin = min(cont_shi_d(:));  
   c_int=cmin:(cmax-cmin)/1000 :cmax;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     contour(x,z,cont_shi_d,10)    
       hold on
              n1=3 ; n2= 50;
       contour(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/2:c_int(n2)])
         n1=10 ; n2= 950;
        contour(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/10:c_int(n2)])
          n1=901 ; n2= 989;
             contour(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/3:c_int(n2)])
%       
                      n1=980 ; n2= 992;
        contour(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])                 
% %           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                 n1=996 ; n2= 998;
      contour(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])
                     n1=994 ; n2= 998;
      contour(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])
                     n1=999 ; n2= 1001;
              contour(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])
% %
 %% %
 colormap(jet)       
   
 return
 function cont(x,z,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(x) ;
   for j=1:N
        cont_shi_d(:,j)=x(j)*s(:,j);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cmax = max(cont_shi_d(:));
 cmin = min(cont_shi_d(:));  
   c_int=cmin:(cmax-cmin)/1000 :cmax;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     contour(x,z,cont_shi_d,10)    
       hold on
              n1=3 ; n2= 50;
       contourf(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/2:c_int(n2)])
         n1=10 ; n2= 950;
        contourf(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/10:c_int(n2)])
          n1=901 ; n2= 989;
             contourf(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/3:c_int(n2)])
%       
                      n1=980 ; n2= 992;
        contourf(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])                 
% %           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                 n1=996 ; n2= 998;
      contourf(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])
                     n1=994 ; n2= 998;
      contourf(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])
                     n1=999 ; n2= 1001;
              contourf(x,z,cont_shi_d,[c_int(n1):(c_int(n2)-c_int(n1))/1:c_int(n2)])
% %

% 
%    colormap(jet)
%           hold on
%           h = colorbar;
%         h.Limits = [min(cont_shi_d(:)) max(cont_shi_d(:))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
 
 
    
%     h2=title(['$\mu= $' num2str(xmu),',   ','$\xi= $' num2str(xxi),',   ',...
%      '$\nu=$' num2str(xnu)]);
%   set(h2,'interpreter','LaTex','Fontsize',12);