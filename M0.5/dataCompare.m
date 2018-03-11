clc
clear
NFFT=8192;

MaX=0,0; MaY=0.0;MaZ=0.0;


hold on
box on 
grid on

px=load('pF.txt');
pF=px(:,1)+1i*px(:,2);
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
tdref=importdata('FDTimePressure1.dat');

%plot(OTime*1000+9.3,-2*f1','r-.','linewidth',1.5);
plot(tdref(:,1)*1000+9,-2*f1','r-.');% ifft from C code
%ref=load('FDTimePressure1.txt');
hold on 

thref=importdata('Subrotatingmonopoletimehistory000.dat');
plot(thref(:,1),thref(:,2),'b-');

axis([170 205 -0.15 0.15])
%legend('C_{ifft}','Ref.');
%


hold on
plot(tdref(:,1)*1000,tdref(:,2),'g-');

legend('ifft_{C}','zhongjie','ref');


Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename1,'-depsc'); 
