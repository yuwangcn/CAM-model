function Enz_v=CAMSVel(t,s)
global ExPara;
global ExVari;
global Ex;
%global Chl_O2;
global O2air;
global CO2air;
global Temp_air0;
%global TempAverage;
global TempDifference;
%global RH;
global I;
global Control_MAL_out;
global Control_MAL_in
global Time_PS;
global Temp_leaf;
 
% global conT
% if conT==1
% Temp_airt=Temp_air0;
% end
% if conT==0
% Temp_airt=Temp_air0+TempDifference*sin((t-3600*4)/43200*pi-pi);
% end
% Temp_leaf=Temp_airt;

Ci=s(1);
Cyto_CO2=s(2);
Cyto_HCO3=s(3);
Cyto_OAA=s(4);
Cyto_PEP=s(5);
Cyto_NADH=s(6);
Cyto_malate=s(7);
Cyto_ATP=s(8);
Vacu_MalicAcid=s(9);
Chl_malate=s(10);
Chl_NADPH=s(11);
Chl_pyruvate=s(12);
Chl_PEP=s(13);
Chl_ATP=s(14);
Chl_CO2=s(15);
Chl_PGA=s(16);
Cyto_C6=s(17);
H2O_EP=s(18);%evaporation
O2i=s(19);
Cyto_O2=s(20);
Chl_O2=s(21);


Vm_1 = 200;% 4.2.1.1
Vm_2o = 0.035;% 4.1.1.31   
Vm_3 = 0.03;% 1.1.1.82       
Vm_4 = 0.03;% 1.1.1.40 
Vm_5 = 0.03;% 2.7.9.1 
Vm_6o = 0.04;% 4.1.1.39
Vm_ATPM = 0.3;
Vm_NADPHM = 0.2; 
Vm_62 =0.03;%mM/s %1.45 E-5;
Vm_MAL_chl=0.08;
Vm_MAL_vin=0.08;
Vm_MAL_vout=0.08;

Jmax=0.3;
gso=0.1;
gm =0.1;%%molm-2 bar-1
Pco2_B=20;
PH_lim=3.5;

if ExPara==1
    gso=gso*ExVari;
end
if ExPara==2
    gm=gm*ExVari;
end
if ExPara==3
    Vm_6o=Vm_6o*ExVari;
end
if ExPara==4
    Vm_2o=Vm_2o*ExVari;
end

if ExPara==5
    PH_lim=PH_lim-1+ExVari;
end

if ExPara==6
    Vm_MAL_vin=Vm_MAL_vin*ExVari;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sc25=3*10^4;%ubarL/mmol ubarL/mmol
So25=81.58;%kpa L /mmol
%Henry's law constant for CO2 and O2
C_CO2=2400;
C_O2=1700;
kH_CO2=1/Sc25*exp(C_CO2*(1/(Temp_leaf+273.15)-1/298));
kH_O2=1/So25*exp(C_O2*(1/(Temp_leaf+273.15)-1/298));
Sc=1/kH_CO2;
So=1/kH_O2;
%Chl_O2a=Chl_O2*kH_O2;%WY2018
%Chl_O2=Chl_O2*kH_O2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.2.1.1        1
KmCO2_1 = 2.8;%2.8;  
Ke_1 =11.2;%20;%;% No unit  k=4.3*10^(-7)  PH=7.5 13.6 11.82
%4.1.1.31       2 
KmHCO3_2 =0.02; KmPEP_2 =0.3;   Kimal_2 =0.3; %WY201802 Kimal_2 =1;KmPEP_2 =1;%WY201810: 1->0.3
f_Km_2o=10;
%1.1.1.82       3
KmNADPH_3 =0.39;  KmOAA_3 =3.23;  KmNADP_3 =0.69;  Kmmal_3 =3.52; %WY201812 KmNADPH_3 =0.024;  KmOAA_3 =0.056;  KmNADP_3 =0.073;  Kmmal_3 =32.0;
Ke_3 =4450.0; % No unit
%1.1.1.40       4 
KmCO2_4 =1.1;  KmNADP_4 =0.0099;  KmNADPH_4 =0.045;  KmPyr_4 =3;  Kmmal_4 =0.35;  Ke_4 =0.051;%KmNADP_4 =0.008 %WY Kmmal0=0.23 Kmmal=0.35//KmNADP_4 =0.0080,KmNADP_4 =0.0099
%2.7.9.1        5 
KiPEP_5 =0.15;  KmATP_5 =0.082;  KmPyr_5 =0.25;
%4.1.1.39       6 
KmCO2_6 =0.0162;  KmO2_6 =0.183;  KmRuBP_6 =0.02;  KiPGA_6 =2.52;  KiFBP_6 =0.04;  KiSBP_6 =0.075;  KiPi_6 =0.9*3;  KiNADPH_6 =0.07*3;%KiNADPH_6 =0.07;KiPi_6 =0.9

Q=0.7;
T0=0.0013;
%3.6.3.14:MChl   ATPM
KmADP_ATPM =0.014;  KmATP_ATPM =0.11;  KmPi_ATPM =0.3;
Ke_ATPM =5.734;        %1/millimolarity
%1.18.1.2:MChl  NADPHM
KmNADP_NADPHM =0.05;  KmNADPH_NADPHM =0.058;  
Ke_NADPHM =502;  % No unit

% Mutase and enolase


KmPGA_62 = 0.08; 
KmPEP_62 = 0.3;
Ke_62=0.4302;%G66 = +0.5; 


%Km_MAL_Cyto=0.1;
Km_MAL_Cyto=2.5;%WY201811
Km_MAl_Vacu=0.1;
Km_MAl_Chl=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Tempreature response of enzymes
R=8.3144598;% m2 kg s-2 K-1 mol-1
KCA25=124*1.03;
Ea_KCA=40.9*1000;
dS_KCA=0.21*1000;
Hd_KCA=64.5*1000;
% % Kh25=0.038;
% % Ea_Kh=95*1000;
Vpmax25=450*1.01;
Ea_Vpmax=65*1000;% WY201810:94.8*1000
% dS_Vpmax=0.25*1000;
% Hd_Vpmax=73.3*1000;
% Kp25=60.5;
% Ea_Kp=27.2*1000;
Vcmax25=0.89;
Ea_Vcmax=78*1000;
% Voc25=0.16*Vcmax25;
% Ea_Voc=55.3*1000;
Kc25=121;
Ea_Kc=64.2*1000;
Ko25=29.2;
Ea_Ko=10.5*1000;
Vm_1 = Vm_1*exp(Ea_KCA*(Temp_leaf+273.15-298.15)/(298.15*R*(Temp_leaf+273.15)))*(1+exp((298.15*dS_KCA-Hd_KCA)/(298.15*R)))/(1+exp(((Temp_leaf+273.15)*dS_KCA-Hd_KCA)/((Temp_leaf+273.15)*R)));
% %Kh(i)= Kh25*exp(Ea_Kh*((Temp_leaf+273.15)-298.15)/298.15/R/(Temp_leaf+273.15));
%Vm_2o= Vm_2o*exp(Ea_Vpmax *((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)))*(1+exp((298.15*dS_Vpmax-Hd_Vpmax)/(298.15*R)))/(1+exp(((Temp_leaf+273.15)*dS_Vpmax-Hd_Vpmax)/((Temp_leaf+273.15)*R)));
Vm_2o= Vm_2o*exp(Ea_Vpmax *((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));%WY201811
Vm_6o= Vm_6o*exp(Ea_Vcmax*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
% Vomax(i)=Voc25*exp(Ea_Voc*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
KmCO2_6= KmCO2_6*exp(Ea_Kc*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
KmO2_6= KmO2_6*exp(Ea_Ko*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));
%%%%%%
%PPDK
Km25=85.4362/1000;
Ea_PPDK=25.998*1000;
KmPyr_5= KmPyr_5*exp(Ea_PPDK*((Temp_leaf+273.15)-298.15)/(298.15*R*(Temp_leaf+273.15)));

Vm_5=Vm_5*(-185.19+83.66*Temp_leaf-1.745*Temp_leaf^2+0.00673*Temp_leaf^3);
%%%MEandMDH
Q10=2;
Vm_3=Vm_3*Q10^((Temp_leaf-25)/10);
Vm_4=Vm_4*Q10^((Temp_leaf-25)/10);


global vstomata;
global vstomata_H2O;

% Tleafandwaterlose = balanceLeafTemperature(0, gs*1000, Temp_air);
% Temp_leaf=Tleafandwaterlose(1);

%H2Oair=6.112*exp(17.62* Temp_leaf/(243.12 + Temp_leaf))*RH;%keep RH constant
%H2Oair=6.112*exp(17.62* 25/(243.12 + 25))*0.6;% keep WVP constant
H2Oair=6.1094*exp(17.625* 10/(243.04 + 10));% keep WVP constant%WY201802
H2Oi=6.1094*exp(17.62* Temp_leaf/(243.04 + Temp_leaf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Chl_NADP=0.5-Chl_NADPH;
Chl_ADP=1.5-Chl_ATP;
Chl_Pi=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod (ceil(t/43200),2)==0%%%12 hour light
%    It= I;%1*sin(t-43200)/43200*pi
     gs=0;%0.001*gs0;%%WY201811
     if t>43200+floor(t/86400)*86400&&t<=44820+floor(t/86400)*86400
     f_Km_2=f_Km_2o*(1-((t-floor(t/86400)*86400)-43200)/1800); 
     Vm_2=Vm_2o;
     %f_Km_2=f_Km_2o;
     %Vm_2=Vm_2o*(1-((t-floor(t/86400)*86400)-43200)/1800); 
    % end
     else%if t>45000+floor(t/86400)*86400
     f_Km_2=0.1*f_Km_2o;
     Vm_2=Vm_2o;
%      Vm_2=0;
     end
     if Ci<=400
     gs=gso;%0.03;    
     %Vm_2=0.2*Vm_2o;
     end
     gs=gs*sqrt(1-(Ci/CO2air)^2)*(H2Oair/H2Oi); %%%WY 201801
     Vm_6=0;
     if t<=68400+floor(t/86400)*86400
     Vm_6=Vm_6o*(0.4+((t-floor(t/86400)*86400)-43200)/42000);
     end
     if t>68400+floor(t/86400)*86400&&t<=79200+floor(t/86400)*86400
     Vm_6= Vm_6o;
     end
     if t>79200+floor(t/86400)*86400%&&t<floor(t/86400)*86400;
     Vm_6=Vm_6o*(1-((t-floor(t/86400)*86400)-79200)/(7200*2));
     end 
     if t==floor(t/86400)*86400;
     Vm_6=Vm_6o*0.5;
     end
     
   Control_MAL_out=0;
   Control_MAL_in=1;  
   if Vm_6>0.5*Vm_6o
    Control_MAL_out=1;
    Control_MAL_in=0;
   end
    if Vacu_MalicAcid<=10;
    Control_MAL_out=0;
    end



else %%%%%dark
    Vm_6=0;
%    It=0;
     gs=gso;
     if t<12600+floor(t/86400)*86400;
     gs=gso*(0.03+0.02*(t-floor(t/86400)*86400)/3600)/0.1;
     end
     gs=gs*sqrt(1-(Ci/CO2air)^2)*(H2Oair/H2Oi); %%%WY 201801

%    Vm_2=Vm_2o*(0.2+0.1*(t-floor(t/86400)*86400)/3600);
    Vm_2=Vm_2o;
    f_Km_2=f_Km_2o*(0.2+0.1*(t-floor(t/86400)*86400)/3600);
    if f_Km_2>f_Km_2o
        f_Km_2=f_Km_2o;
    end
    
   if Vm_2 >Vm_2o
       Vm_2 =Vm_2o;
   end
   if Cyto_C6<10
       Vm_2 =0;
   end
       
   Control_MAL_out=0;
   Control_MAL_in=1;

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if Cyto_C6<100
%    Vm_2=0;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   It=0; 
%   I_0=I;
%   It=I_0*sin(t/43200*pi-pi);
%   if It<0
%       It=0;
%   end
% if t>43200+floor(t/86400)*86400 &&t<=43200+3600+floor(t/86400)*86400;
%    It=I_0*(t/43200*pi-pi); 
% end
% if t>43200+3600+floor(t/86400)*86400&&t<=43200+3600*2+floor(t/86400)*86400;
%    It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*2+floor(t/86400)*86400&&t<=43200+3600*3+floor(t/86400)*86400;
%    It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*3+floor(t/86400)*86400&&t<=43200+3600*4+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*4+floor(t/86400)*86400&&t<=43200+3600*5+floor(t/86400)*86400;
%     It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*5+floor(t/86400)*86400&&t<=43200+3600*6+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*6+floor(t/86400)*86400&&t<=43200+3600*7+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*7+floor(t/86400)*86400&&t<=43200+3600*8+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*8+floor(t/86400)*86400&&t<=43200+3600*9+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*9+floor(t/86400)*86400&&t<=43200+3600*10+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*10+floor(t/86400)*86400&&t<=43200+3600*11+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% if t>43200+3600*11+floor(t/86400)*86400&&t<=43200+3600*12+floor(t/86400)*86400;
%      It=I_0*(t/43200*pi-pi);
% end
% 
% light input

It=0;
if t>43200+floor(t/86400)*86400 &&t<=86400+3600+floor(t/86400)*86400;
   T_0=43200:3600:86400;%43200:3600:86400;%%%%%%%%%%%%%%%
   I_0=I;
   %It = interp1(T_0,I_0,t);
   It = interp1(T_0,I_0,t-floor(t/86400)*86400);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% gs
% Ci
% H2Oi
% xx=(1-(Ci/CO2air)^2)*(H2Oair/H2Oi)
% if t>43200+floor(t/86400)*86400 &&t<=86400+3600+floor(t/86400)*86400;
% if Ci<CO2air
% gs=gs*(1-(Ci/CO2air)^2)*(H2Oair/H2Oi);
% else
% gs=0;
% end
% else
% gs=gs*(1-(Ci/CO2air)^2)*(H2Oair/H2Oi);   
% end


% if gs==0;
%     gs=0.000001;
% end

% % gs=gs*(1-(Ci/CO2air)^2)*(H2Oair/H2Oi);
% % if t>43200+floor(t/86400)*86400 &&t<=86400+3600+floor(t/86400)*86400 && Ci>CO2air
% %     gs=0;
% % end   

vstomata=gs*10^(-3)*(CO2air-Ci);%mmol m-2 s-1bar-1 %gs mol m-2 s-1 0.2*1000/1000000
vstomata_H2O=gs*(H2Oi-H2Oair)*1.6;

global vinf;
vinf=gm*10^(-3)*(Ci-Cyto_CO2*Sc);%%% *Sc added wy201809

%C4 Cycle 5
global v1;
v1=Vm_1*(Cyto_CO2-Cyto_HCO3/Ke_1)/(KmCO2_1+Cyto_CO2);

global v2;
if Cyto_C6>=0&&Cyto_malate<300 %%WY 201707
%v2=Vm_2*Cyto_HCO3*Cyto_PEP/(Cyto_PEP+KmPEP_2*(1+Cyto_malate/(Kimal_2*f_Km_2)))/(Cyto_HCO3+KmHCO3_2);
v2=Vm_2*Cyto_HCO3*Cyto_PEP/(Cyto_PEP+KmPEP_2)/(Cyto_HCO3+KmHCO3_2)*(1-(Cyto_malate/(Cyto_malate+Kimal_2*f_Km_2)));
else
v2=0;
end
%v2=Vm_2*Cyto_HCO3*Cyto_PEP/(Cyto_PEP+KmPEP_2*(1+Cyto_malate/(Kimal_2)))/(Cyto_HCO3+KmHCO3_2);
global v3;
Cyto_NAD=0.2;
if  Cyto_OAA>0.001%%WY 201707
v3=Vm_3*(Cyto_OAA*Cyto_NADH-Cyto_NAD*Cyto_malate/Ke_3)/( KmOAA_3* KmNADPH_3*(1+Cyto_OAA/KmOAA_3+ Cyto_NADH/KmNADPH_3+ Cyto_NAD/KmNADP_3+ Cyto_malate/Kmmal_3+ Cyto_OAA*Cyto_NADH/(KmOAA_3* KmNADPH_3)+ Cyto_NAD*Cyto_malate/(KmNADP_3* Kmmal_3)));
else
v3=0;
end
global v4;

if Chl_malate>=0%&&Chl_NADP>0.0001%WY201707
v4=Vm_4*(Chl_malate*Chl_NADP-Chl_pyruvate*Chl_NADPH*Chl_CO2/Ke_4)/(Kmmal_4*KmNADP_4)/(1+Chl_malate/Kmmal_4+Chl_NADP/KmNADP_4+Chl_pyruvate/KmPyr_4+Chl_NADPH/KmNADPH_4+Chl_CO2/KmCO2_4+Chl_malate*Chl_NADP/(Kmmal_4*KmNADP_4)+Chl_pyruvate*Chl_NADPH/(KmPyr_4*KmNADPH_4)+Chl_pyruvate*Chl_CO2/(KmPyr_4*KmCO2_4)+Chl_NADPH*Chl_CO2/(KmNADPH_4*KmCO2_4)+Chl_pyruvate*Chl_NADPH*Chl_CO2/(KmPyr_4*KmNADPH_4*KmCO2_4));
else
v4=0;
end
if Chl_pyruvate<0 && v4<0 
v4=0;
end
    
global v5;
if Chl_pyruvate>=-0
v5=v4;%Vm_5*Chl_pyruvate*Chl_ATP/(Chl_pyruvate+KmPyr_5*(1+Chl_PEP/KiPEP_5))/(Chl_ATP+KmATP_5);
else
v5=0;
end
% Calvin cycle 12
global v6;
I2=It*0.85*0.85/2;
J=(I2+Jmax-sqrt(((I2+Jmax)^2-4*Q*I2*Jmax)))/(2*Q);
% Ka6_ATP=1;
% if Chl_ATP<=1%WY201802
% Vm_6=Vm_6*(Chl_ATP/Ka6_ATP);
% end
if Chl_ATP>0
Ka6_ATP=0.2;
Vm_6=Vm_6*Chl_ATP/(Chl_ATP+Ka6_ATP);
else
Vm_6=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ac=(Chl_CO2-T0)*Vm_6/(Chl_CO2+KmCO2_6*(1+Chl_O2/KmO2_6));
% Aj=(Chl_CO2-T0)*J/(4*Chl_CO2+8*T0);
% Ao=0;
% v6=min(Ac,Aj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%201811 WY
vomax=0.11*Vm_6;
%Ao=vomax*KmCO2_6/KmO2_6*Chl_O2/Chl_CO2;
T0x=0.5*vomax*KmCO2_6/KmO2_6*Chl_O2;
%Ac=(Chl_CO2-T0x)*Vm_6/(Chl_CO2+KmCO2_6*(1+Chl_O2/KmO2_6));
%Aj=(Chl_CO2-T0x)*J/(4*Chl_CO2+8*T0);
%v6=min(Ac,Aj);
Ac=(Chl_CO2)*Vm_6/(Chl_CO2+KmCO2_6*(1+Chl_O2/KmO2_6));
Ao=(Chl_O2)*vomax/(Chl_O2+KmO2_6*(1+Chl_CO2/KmCO2_6));
v6=Ac;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v6=Ac;

if v6<0
    v6=0;
end

%ATP&NADPH 3
vATP=min(Vm_ATPM,J)*(Chl_ADP*Chl_Pi-Chl_ATP/(Ke_ATPM))/(KmADP_ATPM*KmPi_ATPM*(1+Chl_ADP/KmADP_ATPM+Chl_Pi/KmPi_ATPM+Chl_ATP/KmATP_ATPM+Chl_ADP*Chl_Pi/(KmADP_ATPM*KmPi_ATPM)));
vNADPH=min(Vm_NADPHM,J/2)*(Chl_NADP-Chl_NADPH/Ke_NADPHM)/(KmNADP_NADPHM*(1+Chl_NADP/KmNADP_NADPHM+Chl_NADPH/KmNADPH_NADPHM));



 
%v6=Vm_6*Chl_RuBP*Chl_CO2/((Chl_CO2+KmCO2_6*(1+Chl_O2/KmO2_6))*(Chl_RuBP+KmRuBP_6*(1+Chl_PGA/KiPGA_6+Chl_FBP/KiFBP_6+Chl_SBP/KiSBP_6+Chl_Pi/ KiPi_6)));
%v6=Vm_6*Chl_RuBP*Chl_CO2/((Chl_CO2+KmCO2_6*(1+Chl_O2/KmO2_6))*(Chl_RuBP+KmRuBP_6));

global v62_chl;
v62_chl= Vm_62*Chl_PEP/(KmPEP_62+Chl_PEP);

global vMAL_Vacu;


% if Cyto_malate>0&&Chl_malate>0
vMAL_Chl=Vm_MAL_chl*(Cyto_malate-Chl_malate)/(Km_MAL_Cyto*(1+Cyto_malate/Km_MAL_Cyto+Chl_malate/Km_MAl_Chl));% 
% else
% vMAL_Chl=0;
% end
    

%%%%%Vacule reactions%%%%%
%Ka1=10^(-3.4)
%Ka2=10^(-5.2)
Ka1_MA=10^(-3.4);
Ka2_MA=10^(-5.2);

% Vacu_H=(Ka1_MA*Vacu_MalicAcid)^0.5;
% %Vacu_malate=(Ka1_MA*Ka2_MA*Vacu_MalicAcid)/(Vacu_H*Vacu_H); 
% Vacu_malate=(Ka1_MA*Vacu_MalicAcid)/(Vacu_H*Vacu_H); 
% 
% 
% Ki_H_Vacu=0.316;
% Vm_MAL_Bin=Vm_MAL_B*(1-Vacu_H/Ki_H_Vacu);
% if Vacu_H>0.316%%%%%%-0.9984
%    Control_MAL_in=0;
% end
% %vMAL_Vacu=Vm_MAL_B*(Cyto_malate-Vacu_malate)/(Km_MAL_Cyto*(1+Cyto_malate/Km_MAL_Cyto+Vacu_malate/Km_MAl_Vacu));% 
% vMAL_in=Control_MAL_in*Vm_MAL_Bin*Cyto_malate/(Km_MAL_Cyto+Cyto_malate);
% vMAL_out=Control_MAL_out*Vm_MAL_B*(Vacu_malate-Cyto_malate)/(Km_MAl_Vacu+Vacu_malate);
% 
% % if Vacu_H>0.316;
% % vMAL_in=0;
% % end
%Vacu_H=(Ka1_MA*Vacu_MalicAcid/1000)^0.5*1000;
%Vacu_H=(Ka1_MA*Vacu_MalicAcid)^0.5;
t1=0.60956;
A1=57365.75584;
y0=39.9996;
if Vacu_MalicAcid>40
pH_vacu=-t1*log((Vacu_MalicAcid-y0)/A1);%Luttge&Smith 1984
else
pH_vacu=6;
end
Vacu_H=10^(-pH_vacu)*1000;
%Vacu_malate=((Ka1_MA*Ka2_MA*Vacu_MalicAcid/1000)/(Vacu_H/1000*Vacu_H/1000))*1000; 
%Vacu_malate=(Ka1_MA*Vacu_MalicAcid)/(Vacu_H*Vacu_H); 
t2=1.42528;
A2=-14.23637;
y2=1.23875;
MalateRatio = A2*exp(-pH_vacu/t2) + y2;

if MalateRatio>1
   MalateRatio=1;
end
% if MalateRatio<A2*exp(-3.5/t2) + y2
%     MalateRatio=A2*exp(-3.5/t2) + y2;
% end
if MalateRatio<=0.03
    MalateRatio=0.03;
end

Vacu_malate=Vacu_MalicAcid*MalateRatio;

Ki_H_Vacu=10^(-PH_lim)*1000;%0.316;%pH=3.510^(-3.5)
Vm_MAL_Bin=Vm_MAL_vin*(1-Vacu_H/Ki_H_Vacu);

if Vacu_H>Ki_H_Vacu%%%%%%pH=3.5
   Control_MAL_in=0;
end
%vMAL_Vacu=Vm_MAL_B*(Cyto_malate-Vacu_malate)/(Km_MAL_Cyto*(1+Cyto_malate/Km_MAL_Cyto+Vacu_malate/Km_MAl_Vacu));% 
if Cyto_malate>0&&pH_vacu>=PH_lim
vMAL_in=Control_MAL_in*Vm_MAL_Bin*Cyto_malate/(Km_MAL_Cyto+Cyto_malate);
else
vMAL_in=0;
end

if Vacu_malate>0%&&Cyto_malate>0
vMAL_out=Control_MAL_out*Vm_MAL_vout*(Vacu_malate-Cyto_malate)/(Km_MAl_Vacu+Vacu_malate);
else
vMAL_out=0;
end
if vMAL_out<0
    vMAL_out=0;
end

%vMAL_out=Control_MAL_out*Vm_MAL_out*(Vacu_MalicAcid)/(Km_MAl_Vacu+Vacu_MalicAcid);




vMAL_Vacu=vMAL_in-vMAL_out;


if vMAL_Vacu >0
    vATP_Vacu=2*vMAL_Vacu;
else
    vATP_Vacu=0;
end

 %vleak_Chl=Pco2_B*(Chl_CO2-Cyto_CO2);%%  Bchl_CO2 -> BSC_CO2
 vleak_Chl=Pco2_B*10^(-3)*Sc*(Chl_CO2-Cyto_CO2);%%  Bchl_CO2 -> BSC_CO2 WY201802
 vC6out_PEP=v2;
 vC6out_ATP=0;%Vm_ATPnight*(Chl_ADP*Chl_Pi-Chl_ATP/(Ke_ATPM))/(KmADP_ATPM*KmPi_ATPM*(1+Chl_ADP/KmADP_ATPM+Chl_Pi/KmPi_ATPM+Chl_ATP/KmATP_ATPM+Chl_ADP*Chl_Pi/(KmADP_ATPM*KmPi_ATPM)));

global vrpd; 
if Cyto_C6>=0
vrpd=0.001;%WY201803 0.0006 % WY201809 0
else
vrpd=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WY201811
vstomata_O2=gs*(O2air-O2i)*1.17*10000;
vinfO2=gm*10^(-3)*1.17*(O2i-Cyto_O2*So)*10000;
vO2_chl=Pco2_B*10^(-3)*So*1.17*(Chl_O2-Cyto_O2);%%  Bchl_CO2 -> BSC_CO2 WY201802

% vstomata_O2=0;%gs*(O2air-O2i)*1.17*10000;
% vinfO2=0;%gm*10^(-3)*1.17*(O2i-Cyto_O2*So)*10000;
% vO2_chl=0;%Pco2_B*10^(-3)*So*1.17*(Chl_O2-Cyto_O2);%%  Bchl_CO2 -> BSC_CO2 WY201802
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Enz_v;
Enz_v=zeros(18,1); 
Enz_v(1)=vstomata;
Enz_v(2)=vinf;
Enz_v(3)=v1;
Enz_v(4)=v2;
Enz_v(5)=v3;
Enz_v(6)=v4;
Enz_v(7)=v5;
Enz_v(8)=v6;
Enz_v(9)=v62_chl;
Enz_v(10)=vMAL_Chl;
Enz_v(11)=vMAL_Vacu;
Enz_v(12)=vleak_Chl;
Enz_v(13)=vATP;
Enz_v(14)=vNADPH;
Enz_v(15)=vC6out_PEP;
Enz_v(16)=vC6out_ATP;
Enz_v(17)=vstomata_H2O;
Enz_v(18)=vrpd;
Enz_v(19)=Ao;
Enz_v(20)=vstomata_O2;
Enz_v(21)=vinfO2;
Enz_v(22)=vO2_chl;
Enz_v=real(Enz_v);
if (t > Time_PS)
    Time_PS = Time_PS + 1;
%    PS_OLD_TIME = t;
end
global PS_VEL;
PS_VEL(1,Time_PS)=t;
PS_VEL(2,Time_PS)=vstomata;
PS_VEL(3,Time_PS)=vinf;
PS_VEL(4,Time_PS)=v1;
PS_VEL(5,Time_PS)=v2;
PS_VEL(6,Time_PS)=v3;
PS_VEL(7,Time_PS)=v4;
PS_VEL(8,Time_PS)=v5;
PS_VEL(9,Time_PS)=v6;
PS_VEL(10,Time_PS)=v62_chl;
PS_VEL(11,Time_PS)=vMAL_Chl;
PS_VEL(12,Time_PS)=vMAL_Vacu;
PS_VEL(13,Time_PS)=vleak_Chl;
PS_VEL(14,Time_PS)=vATP;
PS_VEL(15,Time_PS)=vNADPH;
PS_VEL(16,Time_PS)=vC6out_PEP;
PS_VEL(17,Time_PS)=vC6out_ATP;
PS_VEL(18,Time_PS)=vstomata_H2O;
PS_VEL(19,Time_PS)=vrpd;
PS_VEL(20,Time_PS)=Ao;
PS_VEL(21,Time_PS)=vstomata_O2;
PS_VEL(22,Time_PS)=vinfO2;
PS_VEL(23,Time_PS)=vO2_chl;%%  Bchl_CO2 -> BSC_CO2 WY201802
PS_VEL=real(PS_VEL);
end
