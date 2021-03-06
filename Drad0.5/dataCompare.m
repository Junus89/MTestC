clc
clear
close all
clc

format long

%
set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultaxesLineWidth',1)
set(0,'DefaultaxesFontSize',15)
%


NFFT=8192;

MaX=0.0; MaY=0.0;MaZ=0.0;
OmegaR = 2.0*pi*27.0;
TR = 2.0*pi/OmegaR;
fR = 1.0/TR;
OmegaM = 500*2.0*pi; 
fM = OmegaM/(2*pi);

Tint=27*TR;
ODT=Tint/NFFT;
OTime = ODT*(0:NFFT-1);



figure(1)
hold on
grid on
f=load('IFFT.txt'); %IFFT data from C code
f1=f(:,1)*NFFT;


thref=importdata('Subrotatingdipoletimehistory000.dat');
plot(thref(:,1),thref(:,2),'k*');
plot(OTime*1000+0.37,2*f1','r-');

legend('Mao et al','Predicted')

xlabel('Time [ms]')
ylabel('Acoustic Pressure [Pa]')

axis([165 205 -1.5 1.5])

set(gca,'XTick',(165:5:205))
set(gcf, 'PaperPositionMode','Auto')   % Use screen size

Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename1,'-depsc'); 



%
%Comparision of pressure spectrum.
%
refSp=importdata('Subrotatingdipolespectra000.dat');
Pspec=importdata('FDPressureSpectrum.txt');
f=Pspec(:,1);
pF=Pspec(:,2);


k = 0;
FNum = NFFT;

for j = 1:FNum
        
    if mod(f(j)-round(fM),round(fR))==0 %?
    
        k = k+1;
        pFN(k) = pF(j);
        fN(k) = f(j);
        
    end
    
end

FNumN = length(refSp);


figure(2)
grid on
hold on
box on

plot(refSp(:,1),refSp(:,2),'k*')
stem(refSp(:,1),pFN(1:FNumN),'rs')

legend('Mao et al','Predicted')

xlabel('Frequency [Hz]')
ylabel('Acoustic Pressure [Pa]')

axis([0 1500 -0.0001 0.14])

set(gcf, 'PaperPositionMode','Auto')   % Use screen size
Filename2 = ['spectrum_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename2,'-depsc'); 


