%Import data_group2.mat

%background parameters
p=101100;
T=23+273.15;
R=287;
Uinf=20;
rho=p/(R*T);
L=0.3;
mju=1.84e-05;
Rel=rho*Uinf*L/mju

%insert lab data
x=DATA(1,1)*0.001;% in the unit of [meter]
z=DATA(:,2)*0.001;% in the unit of [meter]

p_tot=DATA(:,7:28);% total pressure at each z position in the unit of [Pa]
p_tot=p_tot';
p_stat=DATA(:,29:36);% static pressure at each z position in the unit of [Pa]
p_stat=p_stat';
p_tot_Pr=mean(DATA(:,37));
p_stat_Pr=mean(DATA(:,38));
y_tot=[190 170 150 130 110 90 70 63 55 50 43 35 30 23 15 10 -10 -30 -50 -70 -90 -110]*0.001;% in the unit of [meter]
dz_tot=[0 0 0 0 0 0 0 -10 -10 0 -10 -10 0 -10 -10 0 0 0 0 0 0 0]*0.001;%in the unit of [meter]

%total pressure

%initial meshgrid 
[Z,Y]=meshgrid(z,y_tot);
[Z,DZ]=meshgrid(z,dz_tot);
% Z=Z+DZ;
% contourf(Z,Y,p_tot);

%interpolating meshgrid
zq=-100:1:100;zq=zq*0.001;
yq=190:-1:-110;yq=yq*0.001;
[Zq,Yq]=meshgrid(zq,yq);
P_tot=interp2(Z,Y,p_tot,Zq,Yq,'spline');
DZq=interp2(Z,Y,DZ,Zq,Yq,'spline');
Zq=Zq+DZq;
% 
% figure(1)
% contourf(Zq,Yq,P_tot,20);colorbar;



%static pressure
y_stat=[200 160 120 80 40 0 -40 -100]*0.001;
[Zs,Ys]=meshgrid(z,y_stat);
Zs=Zs+0.01;
% contourf(Zs,Ys,p_stat);
P_stat=interp2(Zs,Ys,p_stat,Zq,Yq,'spline');
% 
% figure(2)
% contourf(Zq,Yq,P_stat);
% 
% contourf(Zq,Yq,P_tot-P_stat,20);

U=sqrt(2*(P_tot-P_stat)/rho);%velocity profile
Uinf=sqrt(2*(p_tot_Pr-p_stat_Pr)/rho);
u=U-Uinf;


figure(1)
contourf(Zq,Yq,u,20);colorbar;
title('Velocity deficit profile')
xlabel('Z');
ylabel('Y');
axis([-0.075 0.075 -0.05 0.15]);



%wrong drag coefficient
W=200*0.001;H=300*0.001;S=W*H/5;dz=0.001;dy=0.001;
I1=(p_tot_Pr-P_tot)/(0.5*rho*Uinf^2);
I2=(u/Uinf).^2;
Cd0=(sum(sum(I1*dz*dy))+sum(sum(I2*dz*dy)))/S

%momentum thickness correction
thetal=0.664*L*Rel^(-1/2);
thetat=0.036*L*Rel^(-1/5);
MLl=rho*Uinf^2*thetal;
DCd_MLl=MLl/(5*S*(0.5*rho*Uinf^2));
Cd_lar=Cd0-DCd_MLl
MLt=rho*Uinf^2*thetat;
DCd_MLt=MLt/(5*S*(0.5*rho*Uinf^2));
Cd_turb=Cd0-DCd_MLt


%experimental correction
DCd_exp=sum(I1(:,1)*dy+I1(:,end)*dy)/2*W-sum(I2(:,1)*dy+I2(:,end)*dy)/2*W;
Cd_exp=Cd0-DCd_exp


%exponential fitting correction


%corrected plot
figure(2)
up1=(u(:,1)+u(:,end))/2;
up2=u(:,1);
[w,h]=size(u);
up1=up1.*ones(w,h);
contourf(Zq,Yq,u-up1,20);colorbar;
title('Modified velocity deficit profile')
xlabel('Z');
ylabel('Y');
axis([-0.075 0.075 -0.05 0.15]);

% figure(3)
% contourf(Zq,Yq,u-up2,20);colorbar;




