function CarbonGainAverage=CAMCarbonCal(Ii)
global Result;
global Metaout;
Metaout=zeros(1,17);
global CO2air;
CO2air=400;%ubar
%global Temp_air;
global conT
%conT=0;%if 1, constant T=25; if 0,sine function
% Temp_air=25;
% Ttx=0:1:24;
% Temp_air=23+5*sin(t/43200*pi-pi);
global RH;
RH=0.643;%relative humidity 0.643->20.3*100Pa
%CI=0.005; % mM =150 ubar %inter cellular CO2 concentration input
global I;
I=Ii/1000;
% %I = 1;%mmol m-2 s-1 =2000umol m-2 s-1 % light intensity input
% %I = [100,300,500,700,900,1100,1100,900,700,500,300,100];
% %I=[100,100,100,100,100,100,100,100,100,100,100,100,100];
% % Tx=0:1:12;
% % I=100*sin(Tx/12*pi);
% I =500/1000;
global totallitght;
totallitght=sum(I);
global vrpd;
vrpd=0.001;% respiration
% global Chl_O2;
% Chl_O2= 0.2646;%O2 concentration in BSC chloroplast
global O2air;
O2air=21.2;%Kpa
days=1;
CarbonGaini=zeros(1,20);
global Ini;
Ini=CAMSIni;%Initial values 
CG=CAMCarbonDirve(days,Ini);
CarbonGain=CG(1);
waterlose=CG(3);
CarbonGaini(1)=CarbonGain;
waterlosei(1)=waterlose;
Metain=Metaout;
%Metain(1,1)=0.5*CO2air;
Metaouti(1,:)=Metaout;
for days=2:4
    CG=CAMCarbonDirve(days,Metain);
    CarbonGain=CG(1);
    waterlose=CG(3);
    Metain=Metaout; 
    %Metain(1,1)=0.5*CO2air;
    CarbonGaini(days)=CarbonGain;
    waterlosei(days)=waterlose;
    if CG(2)==0
        break
    end
    Metaouti(days,:)=Metaout;
end

CarbonGainTotal=0;
waterloseTotal=0;
for days=5:6
    CG=CAMCarbonDirve(days,Metain);
    CarbonGain=CG(1);
    waterlose=CG(3);
%     if CarbonGain>15
%          CarbonGain=0;
%      end    
    CarbonGainTotal=CarbonGainTotal+CarbonGain;
    waterloseTotal=waterloseTotal+waterlose;
    Metain=Metaout; 
    %Metain(1,1)=0.5*CO2air;
    CarbonGaini(days)=CarbonGain;
    waterlosei(days)=waterlose;
     if CG(2)==0
         break
     end
Metaouti(days,:)=Metaout;
end
Metaouti(days,:)=Metaout;
global waterloseAverage;
CarbonGainAverage(1)=CarbonGainTotal/2;%gm-2
waterloseAverage=waterloseTotal/2;
CarbonGainAverage(2)=waterloseAverage;
WUE=CG(1)/CG(3);
WUE2=CarbonGainAverage(1)/waterloseAverage;

global Rt;
Rt(1)=CG(1);
Rt(2)=CG(3);
Rt(3)=WUE;
Rt(4)=WUE2;
Rt=Rt';
%end