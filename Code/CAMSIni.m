function con = CAMSIni()
global CO2air;
global O2air;
global totallitght;
Sc=3*10^4;%ubarL/mmol ubarL/mmol
Ci=0.5*CO2air;
Cyto_CO2=0.5*Ci/Sc;
Cyto_HCO3=0.01;
Cyto_OAA=0.02;
Cyto_PEP=1;
Cyto_NADH=0.3;
Cyto_malate=0.5;
Cyto_ATP=0.5;
Vacu_MalicAcid=10;%WY201810:1->10
Chl_malate=0.5;
Chl_NADPH=0.3;
Chl_Pyruvate=0.2;
Chl_PEP=0.2;
Chl_ATP=0.5;
Chl_CO2=0.5*Cyto_CO2;
Chl_PGA=0.1;
% if totallitght<=1.8
%     Cyto_C6=1000*totallitght;
% else
Cyto_C6=5000;%Sugar, 6 carbon
% end
H2O_EP=0;
O2i=O2air;
Cyto_O2=0.2646;
Chl_O2=0.2646;

con=zeros(1,18);
con(1)=Ci;
con(2)=Cyto_CO2;
con(3)=Cyto_HCO3;
con(4)=Cyto_OAA;
con(5)=Cyto_PEP;
con(6)=Cyto_NADH;
con(7)=Cyto_malate;
con(8)=Cyto_ATP;
con(9)=Vacu_MalicAcid;
con(10)=Chl_malate;
con(11)=Chl_NADPH;
con(12)=Chl_Pyruvate;
con(13)=Chl_PEP;
con(14)=Chl_ATP;
con(15)=Chl_CO2;
con(16)=Chl_PGA;
con(17)=Cyto_C6;
con(18)=H2O_EP;%evaporation
con(19)=O2i;
con(20)=Cyto_O2;
con(21)=Chl_O2;
end
