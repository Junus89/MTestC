clc
clear

NFFT=8192;
%fff=importdata('pF_New.txt');
%fff=load('pF388802.txt');
fff=load('pF.txt');
pFf=fff(:,1)+1i*fff(:,2);
ipFf=ifft(pFf,NFFT)*NFFT;
ireal=2*real(ipFf);

ODT=1/NFFT;
OTime = ODT*(0:NFFT-1);
figure(1)
hold on
box on
grid on

plot(OTime(1:NFFT),ireal,'r-');

