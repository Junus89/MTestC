clear
close all
clc
tic

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultaxesLineWidth',1)
set(0,'DefaultaxesFontSize',15)
%
MaX = 0.0;MaY = 0.0;MaZ = 0.0;
%
BNum = 1;                                                       
R = 1.0;                                         

%

NFFT = 32768;
OmegaR = 509.9988074009404
TR = 2*pi/(BNum*OmegaR);
fR=1/TR;
Tint=25*TR;
ODT=Tint/NFFT;
OTime = ODT*(0:NFFT-1);

Tint=27*TR;
ODT=Tint/NFFT;
OTime = ODT*(0:NFFT-1);

%%%%%%%%%%%%%%%%%%%Comparision of pressure time history.

px=load('pF.txt');
pF=px(:,1)+1i*px(:,2);
%pFf=fliplr(pF);
ObserverSThetaNum=2;
% 
OpMUpHalf = zeros(NFFT/2+1,ObserverSThetaNum-1);
OpMHalfConj = zeros(NFFT/2-1,ObserverSThetaNum-1);
OpMLowerHalf = zeros(NFFT/2-1,ObserverSThetaNum-1);

for j = 1:ObserverSThetaNum-1
    
    for k = 1:NFFT/2+1
    
        OpMUpHalf(k,j) = pF(k,j);
        
    end
    
end

for j = 1:ObserverSThetaNum-1
    
    for k = 2:NFFT/2 

        OpMHalfConj(k,j) = conj(OpMUpHalf(k,j));
    
    end
    
end

for j = 1:ObserverSThetaNum-1
    
    for k = 2:NFFT/2
        
        OpMLowerHalf(k,j) = OpMHalfConj(NFFT/2-k+1,j);
    
    end
    
end

OpMFull = [OpMUpHalf;OpMLowerHalf];
     
OpTM = zeros(NFFT,ObserverSThetaNum-1);
     
for j = 1:ObserverSThetaNum-1
    
    OpTM(:,j) = real(ifft(OpMFull(:,j),NFFT)*NFFT);

end

% 
f=load('IFFT.txt'); %IFFT data from C code
f1=2*f(:,2)*NFFT;
f1=f(:,1)*NFFT;


figure(1)
hold on 
grid on
box on

thref=importdata('Suprotatingmonopoletimehistory000.dat');
refz=load('FDTimePressure1.txt');
ref = load('FDPressureSpectrum.txt');

%plot(OTime*1000+3,2*f1','r-.','linewidth',1.5); %ifft from C code
%plot(ref(:,1)*1000+3.1,2*f1','r-.','linewidth',1.5); %ifft from C code


plot(thref(:,1),thref(:,2),'k*','linewidth',1.5); % plotting result from ref.
plot(refz(:,1)*1000+49.2,OpTM(:,1),'r-','linewidth',1.5);
%plot(OTime*1000-11,2*f1,'r-','linewidth',1.5);

legend('Poletti et al','Predicted')

xlabel('{\itt} [ms]')
ylabel('{\itp''} [Pa]');
axis([170 205 -1.5 2])

set(gca,'XTick',(170:5:205))
set(gcf, 'PaperPositionMode','Auto')   % Use screen size

Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename1,'-depsc'); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comparision of pressure spectrum %%%%%%%%%%%%
FRef=importdata('Suprotatingmonopolespectra000N.dat');
FDSP1=importdata('FDPressureSpectrum.txt');
FDSP = FDSP1(2:end,:);FDSP(end+1,:)=FDSP1(end,:);

FNum = NFFT;
OmegaR = 509.9988074009404;
OmegaM = 3141.592653589793;
fR = OmegaR/2/pi;
fM = OmegaM/2/pi;
f = FDSP(:,1);
pF = FDSP(:,2);


for k = 1:FNum
    
    fps = fM-round(fM/fR)*fR;
    
    if abs(f(k)-fps)< 10^-0
        
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
    
    if abs(-f(k)-fns)< 10^-0
        
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
