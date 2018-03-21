clc
clear


NFFT=32768;
%NFFT = 8192;


MaX=0.0; MaY=0.0; MaZ=0.0;


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
fR=1/TR;
Tint=25*TR;
ODT=Tint/NFFT;
OTime = ODT*(0:NFFT-1);

% figure(1)
% ref=load('FDTimePressure1.txt');
% 
% 
% thref=importdata('Suprotatingmonopoletimehistory000.dat');
% 
% plot(thref(:,1),thref(:,2),'ko','linewidth',1.5); % plotting result from ref.
% hold on
% 
% %plot(OTime*1000+13,-f1','r-','linewidth',1.5); %ifft from C code
% %plot(OTime*1000+3,c,'b-.','linewidth',1.5);
% %plot(ref(:,1)*1000+3.1,2*f1','r-.','linewidth',1.5); %ifft from C code
% 
% %hold on 
% %plot(ref(:,1)*1000,ref(:,2)/2,'b'); % plotting result from zhongjie
% 
% 
% ObserverSThetaNum=2;
% OpMUpHalf = zeros(NFFT/2+1,ObserverSThetaNum-1);
% OpMHalfConj = zeros(NFFT/2-1,ObserverSThetaNum-1);
% OpMLowerHalf = zeros(NFFT/2-1,ObserverSThetaNum-1);
% 
% for j = 1:ObserverSThetaNum-1
%     
%     for k = 1:NFFT/2+1
%     
%         OpMUpHalf(k,j) = pF(k,j);
%         
%     end
%     
% end
% 
% for j = 1:ObserverSThetaNum-1
%     
%     for k = 2:NFFT/2 
% 
%         OpMHalfConj(k,j) = conj(OpMUpHalf(k,j));
%     
%     end
%     
% end
% 
% for j = 1:ObserverSThetaNum-1
%     
%     for k = 2:NFFT/2
%         
%         OpMLowerHalf(k,j) = OpMHalfConj(NFFT/2-k+1,j);
%     
%     end
%     
% end
% 
% OpMFull = [OpMUpHalf;OpMLowerHalf];
%      
% OpTM = zeros(NFFT,ObserverSThetaNum-1);
%      
% for j = 1:ObserverSThetaNum-1
%     
%     OpTM(:,j) = real(ifft(OpMFull(:,j),NFFT)*NFFT);
% 
% end
% 
% hold on
% 
% plot(OTime*1000+49.2,OpTM(:,1),'r-.','linewidth',1.5);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% axis([170 205 -1.5 2])
% 
% %legend('ifft_{matlab}','ifft_{C}','zhongjie','ref');
% legend('Poletti et al','Predicted','matlab');
% 
% Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
% print(Filename1,'-depsc'); 
% 

figure(2)
hold on;
box on;
grid on;

Fspec = importdata('FDPressureSpectrum.txt');
refSpec=importdata('Suprotatingdipolespectra000.dat');
%mspec=importdata('FDSpectrum1.txt');
stem(refSpec(:,1),refSpec(:,2),'k*');
stem(Fspec(:,1),Fspec(:,2),'rs');

legend('Mao et al','Predicted')

xlabel('{\itf} [Hz]')
ylabel('{\itp''} [Pa]');

axis([-2000 2000 0 0.2])

set(gcf, 'PaperPositionMode','Auto')   % Use screen size
Filename2 = ['spectrum_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename2,'-depsc'); 

%% % Spectrum
% 
% for k = 1:FNum
%     
%     fps = fM-round(fM/fR)*fR;
%     
%     if abs(f(k)-fps)< 10^-4
%         
%         kps = k;
%         break;
%         
%     end
%             
% end
%         
% rpnum = round(fR/f(2));
% 
% k1 = 1;
% 
% while (kps<=FNum)
%     
%     pFN1(k1) = pF(kps);
%     fN1(k1) = f(kps);
%     
%     kps = kps+rpnum;
%     k1 = k1+1;
%     
% end
% 
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% 
% for k = 1:FNum
%     
%     fns = fM-(round(fM/fR)+1)*fR;
%     
%     if abs(-f(k)-fns)< 10^-4
%         
%         kns = -k;
%         break;
%         
%     end
%             
% end
% 
% k2 = 1;
% 
% while (kns>=-FNum)
%     
%     pFN2(k2) = pF(-kns);
%     fN2(k2) = f(-kns);
%     
%     kns = kns-rpnum;
%     k2 = k2+1;
%     
% end
% 
% fN = [-fN2(end:-1:1),fN1];
% pFN = [pFN2(end:-1:1),pFN1];
% 
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% 
% FNumN = length(FRef);
% 
% figure(2)
% grid on
% hold on
% box on
% 
% stem(FRef(:,1),FRef(:,2),'k*')
% stem(FRef(:,1),pFN(1287:1287+FNumN-1),'ro')
% 
% legend('Poletti et al','Predicted')
% 
% xlabel('{\itf} [Hz]')
% ylabel('{\itp''} [Pa]');
% 
% axis([-2000 2000 0 0.02])
% 
% set(gcf, 'PaperPositionMode','Auto')   % Use screen size
% Filename2 = ['spectrum_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
% print(Filename2,'-depsc'); 
