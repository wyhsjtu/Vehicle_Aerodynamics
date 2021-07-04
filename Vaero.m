clear all, close all, clc

A=importdata('data_group2.mat');
A=A.DATA;
z= A(:,2);
y=linspace(20,-11,22);
[Y,Z]=meshgrid(y,z);
Patm= A(1,3);
T= A(1,4)+273;
PtotPr= Patm + A(:,end-1);
PstatPr= Patm + A(:,end);
R = 286.9;
rho=Patm/R/T;
Lp=0.3;
Wp=0.4;
Hp=0.0115;
mu=1.84e-05;
Wah=389e-3/4.17;
Hah=288e-3/4.17;
S=Wah*Hah; %scale 1:4.17
SS=40*50*10^(-4);

PPPtot= Patm + A(:,7:28);
YPt=[19 17 15 13 11 9 7 6.3 5.5 5 4.3 3.5 3 2.3 1.5 1 -1 -3 -5 -7 -9 -11];
DzPt=[0 0 0 0 0 0 0 -1 -1 0 -1 -1 0 -1 -1 0 0 0 0 0 0 0];
Ptot=0*PPPtot;

for i=1:22
PPtot(:,i)=interp1(z-DzPt(i),PPPtot(:,i),z);
end
Ptot=interp2(YPt,z,PPtot,y,z);
Ptot(isnan(Ptot))=1.0093e05;


PPPstat= Patm + A(:,29:end-2);
YPs=[ 20 16 12 8 4 0 -4 -10];
DzPs=[1 1 1 1 1 1 1 1];

for j=1:8
    PPstat(:,j)=interp1(z-DzPs(j),PPPstat(:,j),z);
end
Pstat=interp2(YPs,z,PPstat,y,z);
Pstat(isnan(Pstat))=1.0068e+05;

Uinf= sqrt(mean(PtotPr - PstatPr)/(1/2*rho));
Re=rho*Uinf*Lp/mu;

U=sqrt((Ptot-Pstat)/(1/2*rho));
u=U-Uinf;

contourf(Y,Z,u,10); colorbar;
 
dz=abs(z(2)-z(1))*10^(-3);
dy=abs(y(2)-y(1))*10^(-3);
Dz=abs(z(end)-z(1))*10^(-3);
Dy=abs(y(end)-y(1))*10^(-3);


%% momentum loses theoretical

THl=0.664*Lp*Re^(-1/2);
THt=0.036*Lp*Re^(-1/5);

DMl=rho*Uinf*THl*Dz % Momentum loses due to plate boundary layer in the measurement section
DMt=rho*Uinf*THt*Dz

Cpl=DMl/Wp/Hp/(1/2*rho*Uinf^2);  %Plate drag coef theoreticaly
Cpt=DMt/Wp/Hp/(1/2*rho*Uinf^2);

%% Total drag Coef

CC=1/S*(1/(1/2*rho*Uinf^2)*(dz*dy)*sum(sum(mean(PtotPr)-Ptot)) - (dz*dy)*sum(sum(u.^2/Uinf^2)))

%% drag coef correction

DC=1/S*(1/(1/2*rho*Uinf^2)*Dz*dy*sum(mean(PtotPr)-Ptot(1,:)) - Dz*dy*sum(u(1,:).^2/Uinf^2));

%% Corrected drag of ahmed body

C=CC-DC

%% By fitting data

F=(PtotPr(8)-Ptot(8,:))/PtotPr(8);

%% cftool data a*exp(-((y_normalized-b)/c).^2)

a=0.001222 ;
b=-0.04448  ;
c=  0.629   ; 
c1=a;
c2=(1/c)^2;
c3=b;

yy=(y-4.5)/9.586;  %normalization
g=a*exp(-((yy-b)/c).^2);

plot(y,g,y,F)

Cd=1/Hah/(1/2*rho*Uinf^2)*PtotPr(8)*sum(g)*dy
%% NEW PLOTTING

x=[0,25,35];
y=[0.4795, 0.4857, 0.4816];

plot(x,y,'o-');
xlabel('Slant angle (degree)');
ylabel('Drag coefficient');
legend('Corrected experimental method');











