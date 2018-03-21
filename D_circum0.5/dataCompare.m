clc
clear
NFFT=8192;

MaX=0.0; MaY=0.0;MaZ=0.0;


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


% 
% figure(1)
% hold on
% grid on
% %tdref=importdata('FDTimePressure1.dat');
% 
% %plot(OTime*1000,c,'k-.','linewidth',2.1);
% %plot(tdref(:,1)*1000+9,2*f1','r-.');% ifft from C code
% %ref=load('FDTimePressure1.txt');
% 
% 
% thref=importdata('Subrotatingdipoletimehistory000.dat');
% plot(thref(:,1),thref(:,2),'ko','linewidth',1.5);
% plot(OTime*1000+0.7,2*f1','r-','linewidth',1.5);
% % 
% axis([165 205 -1.2 1.2])
% legend('Mao et al','Predicted');
% %
% 
% 
% %hold on
% % plot(tdref(:,1)*1000,tdref(:,2),'g-');
% % 
% % legend('ifft_{C}','zhongjie','ref');
% 
% 
% Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
% print(Filename1,'-depsc'); 


figure(1)
hold on
grid on
box on

refSp=importdata('D05circum_spec.txt');
Pspec=importdata('FDPressureSpectrum.txt');
%spec=importdata('FDPressureSpectrum2omega194401fnum8192.txt');


plot(refSp(:,1),refSp(:,2),'k*','linewidth',1.1);
%stem(Pspec(:,1)-0.93,Pspec(:,2),'r','linewidth',1.3);
stem(Pspec(:,1),Pspec(:,2),'rs','linewidth',1.1);
axis([-0.,1500,0,0.2]);

legend('Mao et al','Predicted')

xlabel('{\itf} [Hz]')
ylabel('{\itp''} [Pa]');



Filename1 = ['pressurespec_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename1,'-depsc'); 

