% A kinetic model of CAM photosynthesis

function CarbonGain=CAMCarbonDirve(days,Metain)
global Ex;
Ex=1;
global Time_PS;
Time_PS=1;
% global Control_MAL_out;
% global Control_MAL_in

%%%the environment variables%%%%%
global TempAverage;
global Temp_air;
Temp_air=TempAverage;
global Metaout;
%%%the environment variables%%%%%
global CO2air;
%CI=0.005; % mM =150 ubar %inter cellular CO2 concentration input
global I;
global vrpd;
% global Chl_O2;

%global time;
Ini=CAMSIni;%Initial values 
time=43200*2;%*15;%48*3600;%Simulation Time
%options = odeset('RelTol',1e-2,'AbsTol',1e-5);
%options = odeset('Events',@MyEventFunction);
%[Tt, d]=ode23tb(@CAMSMB, [0, time], Metain, options);
tic;
try
%[Tt, d]=ode15s(@CAMSMB, [0, time], Metain, options);
[Tt, d]=ode15s(@CAMSMB, [0, time], Metain);
%d=abs(d);
catch  % 若运行时间超过设定最大时间，则视为模型不稳定
        d=zeros(1,21);
        Tt=0;
        disp('Time Out!');
end

d=real(d);
[row,col]=size(d);
%a=43200*2;
if Tt(row)<time
CarbonGain(1)=0;%gm-2
CarbonGain(2)=0;
CarbonGain(3)=0;
end
if Tt(row)>=time
tx=find(Tt>=time);
CarbonGain(1)=(d(tx(1),17)-d(1,17))*6;%mmol %*30/1000%;%gm-2
CarbonGain(3)=(d(tx(1),18)-d(1,18));%%mmol
CarbonGain(2)=1;
end
global Result;
Result=[Tt, d];
Metaout=d(row,:);
end


