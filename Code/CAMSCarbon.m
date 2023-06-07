% A kinetic model of CAM photosynthesis
function CarbonGain = CAMSCarbon(light,temp,Tcon,Para,Vari)
%clear all;
% global Ex;
% Ex=enzymeR;
global ExPara;
global ExVari;
ExPara=Para;
ExVari=Vari;
global Time_PS;
Time_PS=1;
% global Control_MAL_out;
% global Control_MAL_in

%%%the environment variables%%%%%
global CO2air;
CO2air=400;%ubar
global Temp_air0;
global TempDifference;
Temp_air0=temp;
TempDifference=5;
global RH;
RH=0.643;%relative humidity
global conT
conT=Tcon;%if 1, constant T=25; if 0,sine function
%CI=0.005; % mM =150 ubar %inter cellular CO2 concentration input
global I;
I0=[1,1,1,1,1,1,1,1,1,1,1,1,1];
I=light*I0/1000;
% % Tx=0:1:12;
% % I=light*sin(Tx/12*pi)/1000; %mmol m-2 s-1 =2000umol m-2 s-1 % light intensity input
%I = 1;%mmol m-2 s-1 =2000umol m-2 s-1 % light intensity input
%I = [100,300,500,700,900,1100,1100,900,700,500,300,100];
% I0=[100,100,100,100,100,100,100,100,100,100,100,100,100];
% Tx=0:1:12;
% I=100*sin(Tx/12*pi);
%I =100/1000;
global vrpd;
vrpd=0.001;% respiration
global Chl_O2;
Chl_O2= 0.2646;%O2 concentration in BSC chloroplast

Ini=CAMSIni;%Initial values 
time=43200*2*7;%*15;%48*3600;%Simulation Time
options = odeset('RelTol',1e-3,'AbsTol',1e-5);
[Tt, d]=ode15s(@CAMSMB, [0, time], Ini);
%[Tt, d]=ode23tb(@CAMSMB, [0, time], Ini, options);
%d=abs(d);
d=real(d);
Result =[Tt,d]; 
n=size(Tt); 
Ci=d(:,1);
Cyto_CO2=d(:,2);
Cyto_HCO=d(:,3);
Cyto_OAA=d(:,4);
Cyto_PEP=d(:,5);
Cyto_NADH=d(:,6);
Cyto_malate=d(:,7);
Cyto_ATP=d(:,8);
Vacu_MalicAcid=d(:,9);
Chl_malate=d(:,10);
Chl_NADPH=d(:,11);
Chl_Pyruvate=d(:,12);
Chl_PEP=d(:,13);
Chl_ATP=d(:,14);
Chl_CO2=d(:,15);
Chl_PGA=d(:,16);

Vm_1 = 200;% 4.2.1.1
Vm_2 = 0.03;% 4.1.1.31   
Vm_3 = 0.03;% 1.1.1.82       
Vm_4 = 0.03;% 1.1.1.40 
Vm_5 = 0.03;% 2.7.9.1 
Vm_6 = 0.03;% 4.1.1.39
gs=0.1;
Sc=3*10^4;%ubarL/mmol ubarL/mmol
gm = 0.1;%%molm-2 bar-1
%4.2.1.1        1
KmCO2_1 = 2.8;%2.8;  
Ke_1 =11.2;%20;%;% No unit  k=4.3*10^(-7)  PH=7.5 13.6 11.82
%4.1.1.31       2 
KmHCO3_2 =0.02; KmPEP_2 =1;   Kimal_2 =1;  
%1.1.1.82       3
KmNADPH_3 =0.024;  KmOAA_3 =0.056;  KmNADP_3 =0.073;  Kmmal_3 =32.0;  
Ke_3 =4450.0; % No unit
%1.1.1.40       4 
KmCO2_4 =1.1;  KmNADP_4 =0.0080;  KmNADPH_4 =0.045;  KmPyr_4 =3;  Kmmal_4 =0.23;  Ke_4 =0.051;%0.0344 KmNADP_4 =0.008
%2.7.9.1        5 
KiPEP_5 =0.15;  KmATP_5 =0.082;  KmPyr_5 =0.25;
%4.1.1.39       6 
KmCO2_6 =0.0162;  KmO2_6 =0.183;  KmRuBP_6 =0.02;  KiPGA_6 =2.52;  KiFBP_6 =0.04;  KiSBP_6 =0.075;  KiPi_6 =0.9*3;  KiNADPH_6 =0.07*3;%KiNADPH_6 =0.07;KiPi_6 =0.9
Jmax=0.3;
Q=0.7;
T0=0.0013;

% Mutase and enolase

Vm_62 =0.03;%mM/s %1.45 E-5;
KmPGA_62 = 0.08; 
KmPEP_62 = 0.3;
Ke_62=0.4302;%G66 = +0.5; 

Vm_MAL_B=0.15;
Km_MAL_Cyto=0.1;
Km_MAl_Vacu=0.1;
Km_MAl_Chl=0.1;

Pco2_B=0.2;
Ka1_MA=10^(-3.4);
% 
% for i=1:n
%     
% I2=I*0.85*0.85/2;
% J=(I2+Jmax-sqrt((I2+Jmax)^2-4*I2*Jmax))/(2*Q);
% Ac(i)=(Chl_CO2(i)-T0)*Vm_6/(Chl_CO2(i)+KmCO2_6*(1+Chl_O2/KmO2_6));
% Aj(i)=(Chl_CO2(i)-T0)/(4*Chl_CO2(i)+8*T0);
% A(i)=Ac(i);min(Ac(i),Aj(i));
% %Ainf(i)=gs*10^(-3)*(CO2air-Ci(i));
% Ainf(i)=gm*10^(-3)*(Ci(i)-Cyto_CO2(i)*Sc);
% %Hcon(i)=(Ka1_MA*Vacu_MalicAcid(i))^0.5;
% vleak_Chli(i)=Pco2_B*(Chl_CO2(i)-Cyto_CO2(i));
% vMAL_Chli(i)=Vm_MAL_B*(Cyto_malate(i)-Chl_malate(i))/(Km_MAL_Cyto*(1+Cyto_malate(i)/Km_MAL_Cyto+Chl_malate(i)/Km_MAl_Chl));% 
% 
% Vacu_H(i)=sqrt(Ka1_MA*Vacu_MalicAcid(i));
% Vacu_malate(i)=sqrt(Ka1_MA*Vacu_MalicAcid(i));
% 
% 
% vMAL_Vacu(i)=Vm_MAL_B*(Cyto_malate(i)-Vacu_malate(i))/(Km_MAL_Cyto*(1+Cyto_malate(i)/Km_MAL_Cyto+Vacu_malate(i)/Km_MAl_Vacu));%
% end
[row,col]=size(d);
a=time-43200*2;
tx=find(Tt>=a);
CarbonGain=(d(row,17)-d(tx(1),17))*0.08/1000*6*30;%gm-2
% 
% 
% figure;
% plot(Tt,A*1000,'r');
% hold on;
% plot(Tt,Ainf*1000,'k');
% xlabel('t (s)');
% ylabel('A (\mumol m^-^2 s^-^1)');
% ylim([0,20])
% hold off;
% % figure;
% % plot(Tt,Hcon);
% figure;
% plot(Tt,vleak_Chli*1000,'b');
% hold on
% plot(Tt,vMAL_Chli*1000,'k');
% plot(Tt,vMAL_Vacu*1000,'r');
% ylim([-150,150])
% global PS_VEL;
% figure;
% plot(PS_VEL(1,:),PS_VEL(2,:),'b');
% % hold on;
% % plot(PS_VEL(1,:),PS_VEL(3,:),'k');
% ylim([-0.01,0.02])
% % figure;
% % plot(PS_VEL(1,:),PS_VEL(13,:),'b');
% PS_VEL=PS_VEL';
% [row,col]=size(PS_VEL);
% Leng=find(PS_VEL(:,1)>=86400);
% Num=length(Leng);
% At=zeros(Num,2);
% Numn=row-Num+1;
% At(:,1)=PS_VEL(Numn:row,1)-86400;
% At(:,2)=PS_VEL(Numn:row,2);
% % At(:,1)=PS_VEL(1:Num,1);
% % At(:,2)=PS_VEL(1:Num,2);
% global Enz_v;
% figure;
%  subplot(4,4,1); plot(Tt,d(:,1));
%  title('Ci');
%  subplot(4,4,2); plot(Tt,d(:,2));
%  title('Cyto.CO2');
%  subplot(4,4,3); plot(Tt,d(:,3));
%  title('Cyto.HCO3');
%  subplot(4,4,4); plot(Tt,d(:,4));
%  title('Cyto.OAA');
%  subplot(4,4,5); plot(Tt,d(:,5));
%  title('Cyto.PEP');
%  subplot(4,4,6); plot(Tt,d(:,6));
%  title('Cyto.NADH');
%  subplot(4,4,7); plot(Tt,d(:,7));
%  title('Cyto.malate');
%  subplot(4,4,8); plot(Tt,d(:,8));
%  title('Cyto.ATP');
%  subplot(4,4,9); plot(Tt,d(:,9));
%  title('Vacu.MalicAcid');
%  subplot(4,4,10); plot(Tt,d(:,10));
%  title('Chl.malate');
%  subplot(4,4,11); plot(Tt,d(:,11));
%  title('Chl.NADPH');
%  subplot(4,4,12); plot(Tt,d(:,12));
%  title('Chl.pyruvate');
%  subplot(4,4,13); plot(Tt,d(:,13));
%  title('Chl.PEP');
%  subplot(4,4,14); plot(Tt,d(:,14));
%  title('Chl.ATP');
%  subplot(4,4,15); plot(Tt,d(:,15));
%  title('Chl.CO2');
%  subplot(4,4,16); plot(Tt,d(:,17));
%  title('Cyto_C6');
%  figure;
%  plot(Tt,d(:,18));
end