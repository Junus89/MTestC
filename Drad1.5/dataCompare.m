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

% 
% 
% figure(1)
% hold on
% grid on
% f=load('IFFT.txt'); %IFFT data from C code
% f1=f(:,1)*NFFT;
% 
% 
% thref=importdata('Suprotatingdipoletimehistory000.dat');
% plot(thref(:,1),thref(:,2),'ko');
% plot(OTime*1000+0.37,2*f1','r-');
% 
% legend('Mao et al','Predicted')
% 
% xlabel('Time [ms]')
% ylabel('Acoustic Pressure [Pa]')
% 
% axis([165 205 -1.5 1.5])
% 
% set(gca,'XTick',(165:5:205))
% set(gcf, 'PaperPositionMode','Auto')   % Use screen size
% 
% Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
% print(Filename1,'-depsc'); 
% 
% 

%
%Comparision of pressure spectrum.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comparision of pressure spectrum %%%%%%%%%%%%

FRef=importdata('Suprotatingdipolespectra000.dat');
FDSP=importdata('FDPressureSpectrum.txt');
%FDSP=importdata('FDSpectrum1.txt');

% reff = importdata('FDSpectrum1.txt');
% figure(2)
% hold on
% 
% plot(FDSP(:,1),FDSP(:,2),'g-','linewidth',2.0);
% plot(reff(:,1),reff(:,2),'k-.');
% lenged('man','z');
FNum = NFFT;
OmegaR = 509.9988074009404;
OmegaM = 3141.592653589793;
fR = OmegaR/2/pi;
fM = OmegaM/2/pi;
f = FDSP(:,1);
pF = FDSP(:,2);

for k = 1:FNum
    
    fps = fM-round(fM/fR)*fR;
    
    if abs(f(k)-fps)< 10^-3
        
        kps = k;
        break;
        
    end
            
end
        
rpnum = round(fR/f(2));

k1 = 1;

while (kps<=FNum)
    
    pFN1(k1) = pF(kps);
    fN1(k1) = f(kps);
    
    kps = kps+rpnum;
    k1 = k1+1;
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

for k = 1:FNum
    
    fns = fM-(round(fM/fR)+1)*fR;
    
    if abs(-f(k)-fns)< 10^-4
        
        kns = -k;
        break;
        
    end
            
end

k2 = 1;

while (kns>=-FNum)
    
    pFN2(k2) = pF(-kns);
    fN2(k2) = f(-kns);
    
    kns = kns-rpnum;
    k2 = k2+1;
    
end

fN = [-fN2(end:-1:1),fN1];
pFN = [pFN2(end:-1:1),pFN1];

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

FNumN = length(FRef);


figure(2)
grid on
hold on
box on

stem(FRef(:,1),FRef(:,2),'k*')
stem(FRef(:,1),pFN(1287:1287+FNumN-1),'rs')

legend('Poletti et al','Predicted')

xlabel('{\itf} [Hz]')
ylabel('{\itp''} [Pa]');

axis([-2000 2000 0 0.02])

set(gcf, 'PaperPositionMode','Auto')   % Use screen size
Filename2 = ['spectrum_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename2,'-depsc');


