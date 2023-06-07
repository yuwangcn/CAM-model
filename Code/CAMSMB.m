function EnzMBs = CAMSMB(t,s)  
% global vrpd; 
% vrpd=0.001;%
global conT
global Temp_leaf;
global TempDifference;
global Temp_air0;
if conT==1
Temp_airt=Temp_air0;
end
if conT==0
Temp_airt=Temp_air0+TempDifference*sin((t-3600*4)/43200*pi-pi);
end
Temp_leaf=Temp_airt;
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

Enz_v=CAMSVel(t,s);
vstomata=Enz_v(1);
vinf=Enz_v(2);
v1=Enz_v(3);
v2=Enz_v(4);
v3=Enz_v(5);
v4=Enz_v(6);
v5=Enz_v(7);
v6=Enz_v(8);
v62_chl=Enz_v(9);
vMAL_Chl=Enz_v(10);
vMAL_Vacu=Enz_v(11);
vleak_Chl=Enz_v(12);
vATP=Enz_v(13);
vNADPH=Enz_v(14);
vC6out_PEP=Enz_v(15);
vC6out_ATP=Enz_v(16);
vstomata_H2O=Enz_v(17);
vrpd=Enz_v(18);
Ao=Enz_v(19);
vstomata_O2=Enz_v(20);
vinfO2=Enz_v(21);
vO2_chl=Enz_v(22);
% Volume of compartments

%Assume leaf thickness 200um  1m2 =0.2L         2500um 1m2=2L
%mesophyll cell 50% 0.1L airspace 10% 0.02L     mesophyll cell 80% 
                                                %airspace 5%  0.125L

%    vacule 70%          0.07L                   vacule 90%   1.8L
%    cloroplast 20%     0.02L                    cloroplast 5%  0.1L
%    cytosol 9%           0.009L                 cytosol 4% 0.08L
%    per 1%       0.001L                         per 1% 0.02L

VolCyto=0.08;
VolChl=0.1;
Volper=0.02;
VolVacu=1.8;
VolAirSpace=0.125;
% ODEs
%Delta_Ci=(vstomata-vinf)/VolAirSpace*3*10^4;% WY201809
Delta_Ci=(vstomata-vinf)/VolAirSpace*Sc;% WY201811
%Delta_Cyto_CO2=(vinf-v1+vleak_Chl+vrpd)/VolCyto;
Delta_Cyto_CO2=(vinf-v1+vleak_Chl+vrpd+v2/27*6)/VolCyto;
Delta_Cyto_HCO3=(v1 - v2)/VolCyto;
Delta_Cyto_OAA=(v2 - v3)/VolCyto;
Delta_Cyto_PEP= (v2 - v2)/VolCyto;%
Delta_Cyto_NADH= 0;%
Delta_Cyto_malate=(v3 - vMAL_Vacu -vMAL_Chl)/VolCyto;
Delta_Cyto_ATP=0;
Delta_Vacu_MalicAcid=(vMAL_Vacu)/VolVacu;
Delta_Chl_malate=(vMAL_Chl - v4)/VolChl;
Delta_Chl_NADPH=(vNADPH +v4-2*v6)/VolCyto;
Delta_Chl_pyruvate=(v4 - v5)/VolChl;
Delta_Chl_PEP=(v5 -v5)/VolChl;
Delta_Chl_ATP=(vATP-2*v5-3*v6)/VolCyto;
Delta_Chl_CO2=(v4 - v6 - vleak_Chl)/VolChl;
Delta_Chl_PGA=0;%(v6/3+v62_chl)/VolChl;
%Delta_Cyto_C6=((v6/6+v62_chl/2)-(v2/2+v2/27+vrpd/6))/VolCyto;
Delta_Cyto_C6=((v6/6+v5/2)-(v2/2+v2/27+vrpd/6+0.5*Ao/6));
Delta_H2O_EP=vstomata_H2O;
% WY 201803 CO2 release in cytosol
Delta_Cyto_CO2=(vinf-v1+vleak_Chl+vrpd+v2/27*6+v4)/VolCyto;
Delta_Chl_CO2=(- v6 - vleak_Chl)/VolChl;
% WY 201811 Photorespiration
% Delta_O2i=(vstomata_O2-vinfO2)/VolAirSpace*So;% WY201809
% Delta_Cyto_O2=(vinfO2+vO2_chl-vrpd-v2/27*6)/VolCyto;
% Delta_Chl_O2=(-1.5*Ao + v6 -vO2_chl)/VolChl;
Delta_Chl_ATP=(vATP-2*v5-3*v6-Ao)/VolCyto;
Delta_Cyto_CO2=(vinf-v1+vleak_Chl+vrpd+v2/27*6+v4+0.5*Ao)/VolCyto;

Delta_O2i=0;
Delta_Cyto_O2=0;
Delta_Chl_O2=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global EnzMBs; 
EnzMBs=zeros(18,1);

EnzMBs(1)=Delta_Ci;
EnzMBs(2)=Delta_Cyto_CO2;
EnzMBs(3)=Delta_Cyto_HCO3;
EnzMBs(4)=Delta_Cyto_OAA;
EnzMBs(5)=Delta_Cyto_PEP;%
EnzMBs(6)=Delta_Cyto_NADH;%
EnzMBs(7)=Delta_Cyto_malate;
EnzMBs(8)=Delta_Cyto_ATP;%
EnzMBs(9)=Delta_Vacu_MalicAcid;
EnzMBs(10)=Delta_Chl_malate;
EnzMBs(11)=Delta_Chl_NADPH;
EnzMBs(12)=Delta_Chl_pyruvate;
EnzMBs(13)=Delta_Chl_PEP;
EnzMBs(14)=Delta_Chl_ATP;
EnzMBs(15)=Delta_Chl_CO2;
EnzMBs(16)=Delta_Chl_PGA;
EnzMBs(17)=Delta_Cyto_C6;
EnzMBs(18)=Delta_H2O_EP;
EnzMBs(19)=Delta_O2i;
EnzMBs(20)=Delta_Cyto_O2;
EnzMBs(21)=Delta_Chl_O2;



