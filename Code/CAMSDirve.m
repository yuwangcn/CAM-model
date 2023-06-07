% A kinetic model of CAM photosynthesis

clear all;
global ExPara;
global ExVari;
ExPara=0;
ExVari=1;
global Time_PS;
Time_PS=1;
% global Control_MAL_out;
% global Control_MAL_in

%%%the environment variables%%%%%
global CO2air;
CO2air=400;%ubar
% global Chl_O2;
% Chl_O2= 0.2646;%O2 concentration in BSC chloroplast
global O2air;
O2air=21.2;%Kpa
global Temp_air0;
global TempDifference;
Temp_air0=25;% Average temperature 
TempDifference=5;% Temperature difference
global RH;
RH=0.60;%relative humidity
global conT
conT=1;%if 1, constant T=25; if 0,sine function
%CI=0.005; % mM =150 ubar %inter cellular CO2 concentration input
global I;
I0=[1,1,1,1,1,1,1,1,1,1,1,1,1];
I=1160*I0/1000/50*40;%40 mol m-2 day-1
% % Tx=0:1:12;
% % I=1847*sin(Tx/12*pi)/1000/50; %mmol m-2 s-1 =2000umol m-2 s-1 % light intensity input
global vrpd;
vrpd=0.001;% respiration 1umol m-2 s-1


Ini=CAMSIni;%Initial values 
time=43200*2*7;%*15;%48*3600;%Simulation Time
%options = odeset('RelTol',1e-2,'AbsTol',1e-3);
[Tt, d]=ode15s(@CAMSMB, [0, time], Ini);

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
KmCO2_4 =1.1;  KmNADP_4 =0.0080;  KmNADPH_4 =0.045;  KmPyr_4 =3;  Kmmal_4 =0.23;  Ke_4 =0.051;%0.0344 KmNADP_4 =0.008 %WY Kmmal0=0.23 Kmmal=0.35//KmNADP_4 =0.0080,KmNADP_4 =0.0099
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

global PS_VEL;

figure;
plot(PS_VEL(1,:),PS_VEL(2,:),'b.');
hold on;
% plot(PS_VEL(1,:),PS_VEL(3,:),'k');
ylim([-0.01,0.02])
% figure;
% plot(PS_VEL(1,:),PS_VEL(13,:),'b');
PS_VEL=PS_VEL';
[row,col]=size(PS_VEL);
Leng=find(PS_VEL(:,1)>=86400);
Num=length(Leng);
At=zeros(Num,2);
Numn=row-Num+1;
At(:,1)=(PS_VEL(Numn:row,1)-86400)/3600;
At(:,2)=PS_VEL(Numn:row,2)*1000;
At(:,3)=PS_VEL(Numn:row,20)*1000;
At(:,4)=PS_VEL(Numn:row,9)*1000;
At2(:,1)=(Tt-86400)/3600;
At2(:,2)=d(:,9);%concentration *1.8;
figure;
plot(At(:,1),At(:,2));
hold on;
plot(At2(:,1),At2(:,2)/100);
plot(At(:,1),At(:,3),'b.');
plot(At(:,1),At(:,4),'y.');
% plot(At(:,1),At(:,2)+0.5*At(:,3),'r.');
ylim([0,20]);
;
figure;
 subplot(4,5,1); plot(Tt,d(:,1));
 title('Ci');
 subplot(4,5,2); plot(Tt,d(:,2));
 title('Cyto.CO2');
 subplot(4,5,3); plot(Tt,d(:,3));
 title('Cyto.HCO3');
 subplot(4,5,4); plot(Tt,d(:,4));
 title('Cyto.OAA');
 subplot(4,5,5); plot(Tt,d(:,5));
 title('Cyto.PEP');
 subplot(4,5,6); plot(Tt,d(:,6));
 title('Cyto.NADH');
 subplot(4,5,7); plot(Tt,d(:,7));
 title('Cyto.malate');
 subplot(4,5,8); plot(Tt,d(:,8));
 title('Cyto.ATP');
 subplot(4,5,9); plot(Tt,d(:,9));
 title('Vacu.MalicAcid');
 subplot(4,5,10); plot(Tt,d(:,10));
 title('Chl.malate');
 subplot(4,5,11); plot(Tt,d(:,11));
 title('Chl.NADPH');
 subplot(4,5,12); plot(Tt,d(:,12));
 title('Chl.pyruvate');
 subplot(4,5,13); plot(Tt,d(:,13));
 title('Chl.PEP');
 subplot(4,5,14); plot(Tt,d(:,14));
 title('Chl.ATP');
 subplot(4,5,15); plot(Tt,d(:,15));
 title('Chl.CO2');
 subplot(4,5,16); plot(Tt,d(:,17));
 title('Cyto_C6');
  subplot(4,5,17); plot(Tt,d(:,18));
 title('H2O_EP');
  subplot(4,5,18); plot(Tt,d(:,19));
 title('O2i');
  subplot(4,5,19); plot(Tt,d(:,20));
 title('Cyto_O2');
  subplot(4,5,20); plot(Tt,d(:,21));
 title('Chl_O2');
%  
 
%  figure;
%  plot(Tt,d(:,18));
 [row,col]=size(d);
 R1=[0,0];
 R1(1)=(d(row,17)-5000)*30/1000*6/7;
 R1(2)=(d(row,17)-5000)/d(row,18)*6;
 Numt=[0;0;0;0;0;0;0];
  for i=1:row
     if Tt(i)>=43200*2*1
         Numt(1)=i;
         break;
     end
  end
  for i=1:row
     if Tt(i)>=43200*2*2
         Numt(2)=i;
         break;
     end
  end
  for i=1:row
     if Tt(i)>=43200*2*3
         Numt(3)=i;
         break;
     end
  end
  for i=1:row
     if Tt(i)>=43200*2*4
         Numt(4)=i;
         break;
     end
  end
  for i=1:row
     if Tt(i)>=43200*2*5
         Numt(5)=i;
         break;
     end
 end
 for i=1:row
     if Tt(i)>=43200*2*6
         Numt(6)=i;
         break;
     end
 end
  for i=1:row
     if Tt(i)>=43200*2*7
         Numt(7)=i;
         break;
     end
 end
 %Carbongainperday=(d(row,17)-d(Numt(6),17))*6
 Carbongainperday=[0,0;0,0;0,0;0,0;0,0;0,0;0,0];
 for i=1:7
     if i==1
         Carbongainperday(i,1)=(d(Numt(i),17)-d(1,17))*6*30/1000;
         Carbongainperday(i,2)=(d(Numt(i),17)-d(1,17))/(d(Numt(i),18)-d(1,18))*6;
         
     end
     if i>1
     Carbongainperday(i)=(d(Numt(i),17)-d(Numt(i-1),17))*6*30/1000;
     Carbongainperday(i,2)=(d(Numt(i),17)-d(Numt(i-1),17))/(d(Numt(i),18)-d(Numt(i-1),18))*6;
     end
 end
  Carbongainperday

 
 