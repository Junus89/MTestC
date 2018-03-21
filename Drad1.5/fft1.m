clc
clear

NFFT=32768;
%fff=importdata('pF_New.txt');
%fff=load('pF388802.txt');
fff=load('pF.txt');
pFf=fff(:,1)+1i*fff(:,2);
ipFf=ifft(pFf,NFFT)*NFFT;
ireal=2*real(ipFf);

thref=importdata('Suprotatingmonopoletimehistory000.dat');

ref=importdata('FDTimePressure1.txt');
OmegaR=1.5*340;
TR=2*pi/OmegaR;
Tint=25*TR;
ODT=Tint/NFFT;
OTime = ODT*(0:NFFT-1);
figure(1)
hold on
box on
grid on


plot(ref(:,1)*1000,ireal,'r-');
plot(thref(:,1),thref(:,2),'g-');
%xlim([0,0.3]);
plot(ref(:,1)*1000,0.5*ref(:,2),'k-.');
axis([170 205 -1.5 2])
legend('me','zhongjie','ref');



figure(2)

hold on
grid on
box on

pF=importdata('FDPressureSpectrum.txt');
%mref=importdata('Suprotatingmonopolespectra000N.dat');
%mref=importdata('FDSpectrum1.txt');
stem(pF(:,1),pF(:,2),'r');
%stem(mref(:,1),mref(:,1),'k-.');
hold on
pFref=importdata('Suprotatingmonopolespectra000N.dat');
stem(pFref(:,1),pFref(:,2),'k*');
xlim([-2000,2000])
%legend('com','ref');
%%
pFref=importdata('Suprotatingmonopolespectra000N.dat');
stem(pFref(:,1),pFref(:,2),'k*');
hold on
pF=importdata('FDPressureSpectrum.txt');
%mref=importdata('Suprotatingmonopolespectra000N.dat');
%mref=importdata('FDSpectrum1.txt');
stem(pF(:,1),pF(:,2),'r');
%stem(mref(:,1),mref(:,1),'k-.');
xlim([-2000,2000])
%legend('com','ref');
