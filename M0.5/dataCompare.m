clc
clear
NFFT=8192;

MaX =0.0; MaY=0.0;MaZ=0.0;


hold on
box on 
grid on

px=load('pF.txt');
pF=px(1:end-12,1)+1i*px(1:end-12,2);
%pF=px(:,1)+1i*px(:,2);
%pFf=fliplr(pF);

f=load('IFFT.txt'); %IFFT data from C code
%f1=2*f(:,2)*NFFT;
f1=f(:,1)*NFFT;


%c=real(ifft(pF,NFFT)*NFFT);
c=2*real(ifft(pF,NFFT)*NFFT);

OmegaR=0.5*340;
TR=2*pi/OmegaR;
Tint=27*TR;
ODT=Tint/NFFT;
OTime = ODT*(0:NFFT-1);


figure(1)
hold on
box on
grid on

thref=importdata('Subrotatingmonopoletimehistory000.dat');
plot(thref(:,1),thref(:,2),'ko','linewidth',1.5);

%tdref=importdata('FDTimePressure1.dat');

plot(OTime*1000+9.33,-2*f1','r-','linewidth',1.5);
%plot(OTime*1000+9,-c','r-.');% ifft from C code
%ref=load('FDTimePressure1.txt');
 



axis([170 205 -0.15 0.15])
%legend('C_{ifft}','Ref.');
%



legend('Poletti et al','Predicted');


Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename1,'-depsc'); 

figure(2)
hold on
box on
grid on

refS=importdata('Subrotatingmonopolespectra000N.txt');
Fspec=importdata('FDPressureSpectrum.txt');
stem(Fspec(:,1),Fspec(:,2),'r-','linewidth',1.5);
stem(refS(:,1),refS(:,2),'g-.','linewidth',1.5);
axis([0,1500,0,0.02]);
