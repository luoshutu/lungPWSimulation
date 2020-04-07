fd  = sin(2*pi*(0:0.00001:8));
gg = sin(2*pi*1000*(0:0.00001:8));
lp = sin(2*pi*1000*(0:0.00001:8) + fd);

figure
plot(fd,'r')
figure
plot(lp,'b')
hold on
plot(gg,'r')
axis([0 800 -1 1])

ttp = asin(lp) - asin(gg);
figure
plot(ttp)
% hd   = design(fdesign.lowpass('N,F3dB',16,100,100000),'butter');
% ttp_fil = filter(hd,ttp);
% figure
% plot(ttp_fil)

nn  = length(ttp);
t = 0:nn-1;

Data_I = cos(2*pi*1000/100000*t).*(ttp);
Data_Q = sin(2*pi*1000/100000*t).*(ttp);
figure
plot(Data_I,'r')
figure
plot(Data_Q,'b')

% ÂË²¨
hd   = design(fdesign.lowpass('N,F3dB',16,100,100000),'butter');
Data_I_fil = filter(hd,[Data_I zeros(1,8)]);
Data_Q_fil = filter(hd,[Data_Q zeros(1,8)]);

Data_Amp = sqrt((Data_I_fil).*(Data_I_fil) + (Data_Q_fil).*(Data_Q_fil));
% Data_Amp = sqrt((Data_I).*(Data_I) + (Data_Q).*(Data_Q));

figure
plot(Data_Amp)

Data_e = Data_I_fil + i.*Data_Q_fil;
figure
plot(abs(Data_e))
