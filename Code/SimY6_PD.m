
clear all;

Tmean=[12.2,13.9,17.9,22.1,26.9,32.6,34.7,34.0,30.5,23.8,17.1,11.5];
Tamp=[7.5,7.6,8.3,8.3,8.1,8.0,6.5,6.4,7.0,7.5,7.4,6.9];

for i=1:12
Nofile=num2str(i);
str1='PPF_AgaveY6_PD';
str3='Y6_PD';
str4='test';
str2='.txt';
SC1=[str1,Nofile,str2];
SC2=[str3,Nofile,str2];
SC3=[str3,Nofile,str4,str2];
CAMCanopyCal(SC1,SC2,SC3,0,1,Tmean(i),Tamp(i));
end


