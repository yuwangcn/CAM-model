% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% CAMCanopyCal('PPF_AgaveY6_TD1.txt', 'Y6_TD1.txt','Y1_TD1test.txt',0,0,25,0,25);


clear all;

Tmean=[18.4,19.9,22.0,24.3,25.6,24.9,23.2,23.3,22.8,22.3,20.3,18.7];
Tamp=[7.3,7.9,8.4,8.7,8.3,6.3,5.1,5.1,4.8,5.6,6.7,7.2];

for i=1:12
Nofile=num2str(i);
str1='PPF_AgaveY6_TD';
str3='Y6_TD';
str4='test';
str2='.txt';
SC1=[str1,Nofile,str2];
SC2=[str3,Nofile,str2];
SC3=[str3,Nofile,str4,str2];
CAMCanopyCal(SC1,SC2,SC3,0,1,Tmean(i),Tamp(i));
end


