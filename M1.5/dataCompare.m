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
ObserverSThetaNum=2;

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


%plot(OTime,OpTM(:,1),'r-')





















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
% plot(OTime*1000+3,2*f1','r-.','linewidth',1.5); %ifft from C code
% %plot(ref(:,1)*1000+3.1,2*f1','r-.','linewidth',1.5); %ifft from C code
% 
% hold on 
% plot(ref(:,1)*1000,ref(:,2)/2,'b','linewidth',1.5); % plotting result from zhongjie
% plot(ref(:,1)*1000,OpTM(:,1),'m','linewidth',1.5);
% plot(ref(:,1)*1000,c,'k-.','linewidth',1.5);
% thref=importdata('Suprotatingmonopoletimehistory000.dat');
% plot(thref(:,1),thref(:,2),'g-','linewidth',1.5); % plotting result from ref.
% 
% axis([170 205 -1.5 2])
% 
% %legend('ifft_{matlab}','ifft_{C}','zhongjie','ref');
% legend('ifft_{C}','zhongjie','ifft_{matlab_{complex}}','ifft_{matlab_{simple}}','ref');
% 
% Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
% print(Filename1,'-depsc'); 

figure(2)%compare the spectrum of zhongjie with mine
hold on
box on
grid on
Fspec=importdata('FDPressureSpectrum_w2T194401F32678.txt');
FspecZ=importdata('FDSpectrum1.txt');
stem(Fspec(:,1),Fspec(:,2),'r-','linewidth',1.5);
stem(FspecZ(:,1),FspecZ(:,2),'k-.','linewidth',1.5);
%xlim([0,2000]);
legend('me','zhongjie');

figure(3)
hold on
box on
grid on

%comparing pF and pFT
pFT=importdata('pFT.txt');
plot(OTime,imag(pFT),'r','linewidth',1.6);
plot(OTime,imag(pF),'b-.','linewidth',1.6);
xlim([0,0.02]);
legend('pFT','pF');






% % Spectrum
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
