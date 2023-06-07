function A4totalc = CAMCanopyCal(inputfile, outputfile1,outputfile2,Para,Vari,TempAver,TempDiff)

global ExPara;
global ExVari;
global Temp_air0;
global TempDifference;
ExPara=Para;
ExVari=Vari;
Temp_air0=TempAver;
TempDifference=TempDiff;
global conT
conT=0;%if 1, constant T=25; if 0,sine function
PPFSC = importdata(inputfile);
PPFSC_data =PPFSC.data;
%PPFSC_data = importdata(inputfile);
[row,col]=size(PPFSC_data);
row=row-12;

global Atotal1;
global Atotal2;
global Wtotal1;
global Wtotal2;
global It1;
global It2;
global FA1;
global FA2;
global FW1;
global FW2;
global LA;
%global A4totalc;

XA=100;
XB=300;
YA=100;
YB=300;
%Area=2.56;
Area=4;

PPFDfile1=zeros(row,36);
PPFDfile2=zeros(row,36);
PPFDfile1(:,1:10)=PPFSC_data(1:row,1:10);
PPFDfile1(:,11)=PPFSC_data(1:row,17);
PPFDfile2(:,1:11)=PPFSC_data(1:row,1:11);

for k=1:13
    for m=1:3
    PPFDfile1(:,11+k)=PPFDfile1(:,11+k)+PPFSC_data(1:row,11+m+6*k);%get light intensity
    end
    for n=1:3
    PPFDfile2(:,11+k)=PPFDfile2(:,11+k)+PPFSC_data(1:row,14+n+6*k);%get light intensity
    end
end

    PPFD1(:,1:13)=PPFDfile1(:,12:24);%side1
    PPFD2(:,1:13)=PPFDfile2(:,12:24);%side2
    
FA1=zeros(row,1);
FA2=zeros(row,1);
FW1=zeros(row,1);
FW2=zeros(row,1);
It1=zeros(row,1);
It2=zeros(row,1);
Atotal1=0;
Atotal2=0;
Wtotal1=0;
Wtotal2=0;
X1=PPFDfile1(:,1);
Y1=PPFDfile1(:,2);
for i=1:row  
  %  i
 if X1(i)<=XA||X1(i)>=XB
 FA1(i,1)=0;
 FA2(i,1)=0;
 FW1(i,1)=0;
 FW2(i,1)=0;
 end
 if Y1(i)<=YA||Y1(i)>=YB
 FA1(i,1)=0;
 FA2(i,1)=0;
 FW1(i,1)=0;
 FW2(i,1)=0;
 end
     
if X1(i)>XA&&X1(i)<XB&&Y1(i)>YA&&Y1(i)<YB
i    
Ii1=PPFD1(i,:);%%%%12 hour%%
It1(i)=sum(Ii1);
if It1(i)>460 %100
F1=CAMCarbonCal(Ii1);%CarbonGainAverage
if F1(1)<-86.4%69.1%86.4
    F1(1)=-86.4;
    F1(2)=0;
end
FA1(i,1)=F1(1);
FW1(i,1)=F1(2);
else
%%%%%% Avoid the ODE problem at low light 
FA1(i,1)=0.188*It1(i)-86.4;%FA1(i,1)=-86.4;
if FA1(i,1)<-86.4
    FA1(i,1)=-86.4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FW1(i,1)=0;
end

Ii2=PPFD2(i,:);%%%%12 hour%%
It2(i)=sum(Ii2);
if It2(i)>460%100
F2=CAMCarbonCal(Ii2);%CarbonGainAverage
if F2(1)<-86.4%69.1%86.4
    F2(1)=-86.4;
    F2(2)=0;
end
FA2(i,1)=F2(1);
FW2(i,1)=F2(2);
else
FA2(i,1)=0.188*It2(i)-86.4;%-86.4;
if FA2(i,1)<-86.4
    FA2(i,1)=-86.4;
end
FW2(i,1)=0;
end
end

A1(i)=PPFDfile1(i,11)*FA1(i,1);%%%cm2*g m-2
A2(i)=PPFDfile1(i,11)*FA2(i,1);%%%cm2*g m-2
W1(i)=PPFDfile1(i,11)*FW1(i,1);%%%cm2*g m-2
W2(i)=PPFDfile1(i,11)*FW2(i,1);%%%cm2*g m-2

Atotal1=Atotal1+A1(i); 
Atotal2=Atotal2+A2(i); 
Wtotal1=Wtotal1+W1(i); 
Wtotal2=Wtotal2+W2(i); 

end

LA =zeros(row,1);
global LAI;
LAI=0;
for i=1:row

    if X1(i)>XA&&X1(i)<XB&&Y1(i)>YA&&Y1(i)<YB   
       LA(i)=PPFDfile1(i,11);%%%%%change%%%%
    end    
    if X1(i)<=XA||X1(i)>=XB
       LA(i)=0;
    end
    if Y1(i)<=YA||Y1(i)>=YB
       LA(i)=0;  
    end
    if PPFDfile1(i,11)>25
       LA(i)=0;
     end
     LAI=LAI+LA(i);
end
     


A4totalc=zeros(1,6);
A4totalc(1)=Atotal1/Area/10000;
A4totalc(2)=Atotal2/Area/10000;
A4totalc(3)=Wtotal1/Area/10000;
A4totalc(4)=Wtotal2/Area/10000;
A4totalc(5)=(Atotal1+Atotal2)/(Wtotal1+Wtotal2);
A4totalc(6)=LAI/Area/10000;
Atest=zeros(row,7);
Atest(:,1)=It1;
Atest(:,2)=It2;
Atest(:,3)=FA1;
Atest(:,4)=FA2;
Atest(:,5)=FW1;
Atest(:,6)=FW2;
Atest(:,7)=LA;
dlmwrite(outputfile1,A4totalc,'delimiter','\t','precision', '%.3f');
dlmwrite(outputfile2,Atest,'delimiter','\t','precision', '%.3f');
end

