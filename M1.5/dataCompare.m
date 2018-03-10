clc
clear
NFFT=32768;

hold on
box on 
grid on

px=load('pF.txt');
pF=px(:,1)+1i*px(:,2);
%pFf=fliplr(pF);

f=load('IFFT.txt'); %IFFT data from C code
%f1=2*f(:,2)*NFFT;
f1=f(:,1)*NFFT;
f2=-f(:,1)*NFFT;
f3=fliplr(f2);
f4=fliplr(f1);

%c=real(ifft(pF,NFFT)*NFFT);
c=2*real(ifft(pF,NFFT)*NFFT);

OmegaR=1.5*340;
TR=2*pi/OmegaR;
Tint=25*TR;
ODT=Tint/NFFT;
OTime = ODT*(0:NFFT-1);

figure(1)
ref=load('FDTimePressure1.txt');

%plot(OTime*1000+3,2*f1','r-.','linewidth',1.5); %ifft from C code
plot(ref(:,1)*1000+3.1,2*f1','r-.','linewidth',1.5); %ifft from C code

hold on 
plot(ref(:,1)*1000,ref(:,2)/2,'b'); % plotting result from zhongjie

thref=importdata('Suprotatingmonopoletimehistory000.dat');
plot(thref(:,1),thref(:,2),'g-'); % plotting result from ref.

axis([170 205 -1.5 2])

%legend('ifft_{matlab}','ifft_{C}','zhongjie','ref');
legend('ifft_{C}','zhongjie','ref');
