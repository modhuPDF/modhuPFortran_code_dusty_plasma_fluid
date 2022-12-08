clear all
clc 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Vs=load('Source_v1b_r_15_x_305cy_001_s1_cart.dat');  Vs=Vs(5,:) ;
     omgs=load('Source_w1b_r_15_x_305cy_001_s1_cart.dat');  Ws=omgs(5,:);
          x=load('Xarray1b_r_15_x_305cy_001_cart.dat'); N =length(x); r=x;
          y=load('Yarray1b_r_15_x_305cy_001_cart.dat'); M =length(y); z=y;
  
   w0=load('Output_omg1_mu_1e4_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); 
   p0=load('Output_par1_mu_6e6_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); 
   s0=load('Output_psi1_mu_6e6_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); 
   
%   w0=load('Output_omg1_mu_1e3_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); 
%   w1=load('Output_omg1_mu_1e4_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); %%
%   w2=load('Output_omg1_mu_5e5_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat');%% 
%   w3=load('Output_omg1_mu_2e5_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); %%
%   w4=load('Output_omg1_mu_9e6_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); 
%   w5=load('Output_omg1_mu_6e6_xi_4_nu_3_r_15_x_305cy_toll5_s1_cart.dat'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
 s=s0;
 w=w0;
 p=p0;
xmu=p(3);


figure(1)
clf
%   subplot('Position',[0.2,0.2,0.2,0.5])
cont_shi_d =-s; 
plot_strmlines
% contour(x,z,cont_shi_d,30);  
grid on
    h1=ylabel('$y$');
set(h1,'interpreter','LaTex','Fontsize',15,...
    'FontWeight','bold');
h1=xlabel('$x$');
set(h1,'interpreter','LaTex','Fontsize',15,...
    'FontWeight','bold');
set(gca,'FontSize',13, 'FontName','Arial',...
    'FontWeight','light','Color','w',...
    'LineWidth',1);
 set(gca, 'box', 'on') ;
  h2=title(['$\mu= $' num2str(xmu)]);
 set(h2,'interpreter','LaTex','Fontsize',13);
%  set(gcf,'InnerPosition',[20,20,300,400])
%   set(gcf,'OuterPosition',[20,20,300,400])
 set(gcf,'Position',[20,20,300,400])

 print(gcf,'streamline_6.eps','-depsc2','-r300');

 %%%


return
figure(2)
clf
%    subplot('Position',[0.2,0.2,0.2,0.5])
% cont_shi_d =-w5; plot_strmlines
surf(x,y,-w,'FaceColor','interp',...
'EdgeColor','none',...
'FaceLighting','phong')
% smart_clean('NoWhiteSpace');
% set(gca,'LooseInset',get(gca,'TightInset'))
% view(-30,45)
ylim([y(1),y(M)])
aa= min(w(:));
bb= max(w(:));
% zlim([(aa) (bb)])
% zlim([(aa-0.2*aa) (bb+0.8*bb)])
% zlim([(aa+2.5*aa) (bb+0.2*bb)])
colormap(jet) 
%  grid on
h1=xlabel('$x$');
set(h1,'interpreter','LaTex','Fontsize',15,...
    'FontWeight','bold');
h1=ylabel('$y$');
set(h1,'interpreter','LaTex','Fontsize',15,...
    'FontWeight','bold');
h1=zlabel('$\omega$');
set(h1,'interpreter','LaTex','Fontsize',15,...
    'FontWeight','bold');
grid on
% cax=axis;
% h=text(cax(1)+0.2*(cax(2)-cax(1)),...
%     cax(3)+0.7*(cax(4)-cax(3)),'(a)');
% set(h,'interpreter','LaTex','Fontsize',13);
set(gca,'FontSize',13, 'FontName','Arial',...
    'FontWeight','light','Color','w',...
    'LineWidth',1);
%  set(gcf,'OuterPosition',[20,20,300,400],...
%      'PaperUnits','normalized')
 set(gcf,'InnerPosition',[20,20,300,400])
 
% ss=get(gcf,'position');
% set(gcf,'PaperPosition'[ss]);
%  print(gcf,'vorticity_0','-dpng');
 print(gcf,'vorticity_1','-depsc2','-r300');
%  print(gcf,'vorticity_2','-dpdf');
 