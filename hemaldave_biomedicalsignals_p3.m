clc;
close all;
clear all;
xs=menu('Select the Patient','1) Patient 1-953','2) Patient 2-954','3) Patient 3*955');
switch xs
    
    case 1
        clc;
close all;
clear all;

%% Data load with dlm read function
Data=dlmread('Harding.dat');
p1=Data(2:end,1);
y11=1/600:1/600:size(p1)/600;
Fs=600;
%% performing fft to convert the time domain signal into frequency domain signal.
nfft=length(p1);
nfft2=2.^nextpow2(nfft);
fy=fft(p1,nfft2);
fy=fy(1:nfft2/2);
xfft=Fs.*(0:nfft2/2-1)/nfft2;
%% filter
% to cutoff the power interference signal
f=600; %sampling frequency [Hz]% 
f_cutoff1=42;% lower boundry cutoff frequency
f_cutoff2=42.1;% upper boundry cut off frequenct
order=1;
fnorm1 = f_cutoff1/(f/2);
fnorm2 = f_cutoff2/(f/2);
[b1,a1] = butter(order,[fnorm1,fnorm2],'stop');
 y= filtfilt(b1,a1,p1);% filter output
 % Low pass filter to remove the noise above the 43Hz signal
  f_cutoff1=43;
f_cutoff2=0.2;% cutoff frequency
order=2;
fnorm1 = f_cutoff1/(f/2);
fnorm2 = f_cutoff2/(f/2);
[b1,a1] = butter(order,fnorm1,'low');
 y1= filtfilt(b1,a1,y);% filter output
 %% fft to find the data after filtering
 nfft=length(p1);
nfft2=2.^nextpow2(nfft);
fy1=fft(y,nfft2);
fy1=fy1(1:nfft2/2);
xfft=Fs.*(0:nfft2/2-1)/nfft2;

nfft=length(p1);
nfft2=2.^nextpow2(nfft);
fy2=fft(y1,nfft2);
fy2=fy1(1:nfft2/2);
xfft=Fs.*(0:nfft2/2-1)/nfft2;

%%
select=y1(1:length(y1));% select the section of the data in which we have to perfome peak detection.
[ygh,x]=findpeaks(select,600,'Minpeakdistance',0.75);% 0.75 s window is used 
% to find the Individual RR interval 
i=1:1:size(x)+1;


for i=1:1:size(x)-1
    xi(i)=x(i+1)-x(i);
    hr(i)=60/xi(i);
end
% finding statistic for HRV.
meanNNInterval=sum(xi)/(length(xi)*1000)
heartrate=60/(meanNNInterval*1000)
MinRange=min(xi)
MaxRange=max(xi)
DeviationinNN=std(xi)*1000
% plotting the RR peak time difference to see if the data is uniformly
% sampled or not
k=1:1:length(xi)+1;

for k=1:1:length(xi)-1
    xii(k)=xi(k+1)-xi(k);
end

%% performing frequency domain task
rsf=20;
% as the sampling frequency i am going to resampled the data
Xcs=resample(xi,rsf,1);

% after resampling the data see the change between 2 RR peaks to check the
% uniformity .
k=1:1:length(Xcs)+1;

for k=1:1:length(Xcs)-1
    xiii(k)=Xcs(k+1)-Xcs(k);
end
figure()
% finding the power density spectrum
Fs=20;
% plot psd
N = length(Xcs);
xdft = fft(Xcs);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(Xcs):Fs/2;

% PSD curve analysis

yPower=10*log10(psdx);
% find the given frequency rabge powe by using trapz function.
ulf=abs(trapz(freq(2:10),yPower(2:10)))
lf=abs(trapz(freq(10:35),yPower(10:35)))
hf=abs(trapz(freq(35:92),yPower(35:92)))
tp=ulf+lf+hf
a=lf/hf
% find the max peak in the give range and and find the power by trapz
% function
max=abs(min(yPower(10:35)))
totalpower=abs(trapz(freq(19:21),yPower(19:21)))
coherenceratio=totalpower/(tp-totalpower)
% Fs=30;
% % power density curve
% N = length(xcs);
% xdft = fft(xcs);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/length(xcs):Fs/2;
% figure()
% plot(freq,10*log10(psdx));
% nfft=length(xcs);
% Fs=1;
% % nfft2=2.^nextpow2(nfft);
% fy4=fft(xcs)/nfft;
% fdt=fftshift(fy4);
% fre=(-nfft/2:nfft/2)*rsf/nfft;
% fre=fre(1:end);
% length(fre(end/2:end-1))
% length(abs(fdt(end/2:end)))
% % real plot chris
% figure()
%  plot(fre(end/2:end-1),abs(fdt(end/2:end)));
%% plot
% plot of ECG data before and after resampling.
subplot(3,1,1)
plot(y11,p1)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Patient 1 Data Plot Vs Time');
grid on;
% subplot(3,2,2)
% plot(xfft,abs(fy/max(fy)));
% xlabel('Frequency (in Hz)');
% ylabel('Data(f)');
% grid on;
% title('Single sided FFT plot for Original Data');
subplot(3,1,2)
plot(y11,y)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Filtered data: Bandstop Filter [42-42.1]Hz ');
grid on;
% subplot(3,2,4)
% plot(xfft,abs(fy1/max(fy1)));
subplot(3,1,3)
plot(y11,y1)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Filtered data: low pass filter at 43Hz');
grid on;
% plot the RR peak data.
figure()
findpeaks(select,600,'Minpeakdistance',0.75);
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Finding the peaks');
% plot the HR vs TIme
figure()
plot(x(2:end),hr)
ylabel('Heart Rate');
xlabel(' Time (in seconds)');
title('HR vs Time');
figure()
% plot the sampled RR interval and resample RR interval
subplot(2,1,1)
plot(x(2:end),xi)
xlabel('Time (in seconds)');
ylabel('R-R interval(x(i+1)-x(i))');
title('RR Interval plot');
subplot(2,1,2)
plot(1:1:length(Xcs),Xcs)
xlabel('Time (in seconds)');
ylabel('R-R interval Resampled plot at RSF=20Hz');
title('RR Interval resampled plot AT RSF=20Hz');
figure()
% plot the change between 2 RR interval
subplot(2,1,1)
plot(x(2:end-1),xii);
xlabel('Time (in seconds)');
ylabel('R-R interval Change');
title('RR Interval Change plot');
subplot(2,1,2)
plot(1:1:length(Xcs)-1,xiii);
xlabel('Time (in seconds)');
ylabel('R-R interval Change after Resampling & RSF=20Hz');
title('RR Interval Change  plot after resampling at 20Hz');
figure();
% plotting the PSD
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)');



    
    case 2
        clc;
close all;
clear all;
%% Data load
Data=dlmread('Harding.dat');
p2=Data(2:end,2);
y11=1/600:1/600:size(p2)/600;
Fs=600;
%% performing fft of the original data
nfft=length(p2);
nfft2=2.^nextpow2(nfft);
fy=fft(p2,nfft2);
fy=fy(1:nfft2/2);
xfft=Fs.*(0:nfft2/2-1)/nfft2;
%% filtering
% using band stop filter to remove power interference.
 f=600; %sampling frequency [Hz]% 
f_cutoff1=42;
f_cutoff2=42.15;
order=1;
fnorm1 = f_cutoff1/(f/2);
fnorm2 = f_cutoff2/(f/2);
[b1,a1] = butter(order,[fnorm1,fnorm2],'stop');
 y= filtfilt(b1,a1,p2);
 % used the low pass filter to remove high frequency data which is not
 % needed
  f_cutoff1=43;
f_cutoff2=0.2;
order=2;
fnorm1 = f_cutoff1/(f/2);
fnorm2 = f_cutoff2/(f/2);
[b1,a1] = butter(order,fnorm1,'low');
 y1= filtfilt(b1,a1,y);
 
 %% time domain analysis
 % select the range of the filtred data in which peak detection analysis
 % will happen
 select=y1(1:length(y1));
 % findg the peak
[ygh,x]=findpeaks(select,600,'Minpeakdistance',0.5);
% finding the RR interval and HR
i=1:1:size(x)+1;


for i=1:1:size(x)-1
    xi(i)=x(i+1)-x(i);
    hr(i)=60/xi(i);
end
% performing statistical analysis fo HRV
meanNNInterval=sum(xi)/(length(xi)*1000)
heartrate=60/(meanNNInterval*1000)
MinRange=min(xi)
MaxRange=max(xi)
DeviationinNN=std(xi)*1000

%% Frequency Domain Task
% finding the change in RR interval
k=1:1:length(xi)+1;

for k=1:1:length(xi)-1
    xii(k)=xi(k+1)-xi(k);
end

% after seeing the non unifrom signal , resampling that signal
rsf=20;
Xcs=resample(xi,rsf,1);

% performing this to see the change in RR interval and comparing it with
% without reampled data
k=1:1:length(Xcs)+1;

for k=1:1:length(Xcs)-1
    xiii(k)=Xcs(k+1)-Xcs(k);
end
% performing the PSD
Fs=20;
% plot psd
N = length(Xcs); 
xdft = fft(Xcs);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) *abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(Xcs):Fs/2;

% performing given Frequency power findings

yPower=10*log10(psdx);
ulf=abs(trapz(freq(2:12),yPower(2:12)))
lf=abs(trapz(freq(12:42),yPower(12:42)))
hf=abs(trapz(freq(42:110),yPower(42:110)))
tp=ulf+lf+hf
a=lf/hf
max=abs(min(yPower(12:71)))
totalpower=abs(trapz(freq(47:49),yPower(47:49)))
coherenceratio=totalpower/(tp-totalpower)
%% plotting
% plot the ECG orginal data and filtered data
subplot(3,1,1)

plot(y11,p2)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Patient 1 Data Plot Vs Time');
grid on;
% % subplot(3,2,2)
% plot(xfft,abs(fy/max(fy)));
% xlabel('Frequency (in Hz)');
% ylabel('Data(f)');
% grid on;
% title('Single sided FFT plot for Original Data');
subplot(3,1,2)
plot(y11,y)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Filtered data: Bandstop Filter [42-42.15]Hz ');
grid on;
% subplot(3,2,4)
% plot(xfft,abs(fy1/max(fy1)));
subplot(3,1,3)
plot(y11,y1)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Filtered data: low pass filter at 43Hz');
grid on;
% plot the peak detecttion
figure()
findpeaks(select,600,'Minpeakdistance',0.5);
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Finding the peaks');
grid on
figure()
% plot the HR vs time
plot(x(2:end),hr)
ylabel('Heart Rate');
xlabel(' Time (in seconds)');
title('HR vs Time');
figure()
% plot the RR Interval without resampling and with resampling
subplot(2,1,1)
plot(x(2:end),xi)
xlabel('Time (in seconds)');
ylabel('R-R interval(x(i+1)-x(i))');
title('RR Interval plot');
subplot(2,1,2)
plot(1:1:length(Xcs),Xcs)
xlabel('Time (in seconds)');
ylabel('R-R interval Resampled plot at RSF=20Hz');
title('RR Interval resampled plot AT RSF=20Hz');
figure()
% plot the change in RR interval vs time
subplot(2,1,1)
plot(x(2:end-1),xii);
xlabel('Time (in seconds)');
ylabel('R-R interval Change');
title('RR Interval Change plot');
subplot(2,1,2)
plot(1:1:length(Xcs)-1,xiii);
xlabel('Time (in seconds)');
ylabel('R-R interval Change after Resampling & RSF=20Hz');
title('RR Interval Change  plot after resampling at 20Hz');
figure();
% plot the PSD
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)');
    case 3
%% Data load
Data=dlmread('Harding.dat');
p3=Data(2:end,3);
y11=1/600:1/600:size(p3)/600;
Fs=600;
%% performing fft to analyse the frequecy domain data
nfft=length(p3);
nfft2=2.^nextpow2(nfft);
fy=fft(p3,nfft2);
fy=fy(1:nfft2/2);
xfft=Fs.*(0:nfft2/2-1)/nfft2;


plot(xfft,abs(fy/max(fy)));
xlabel('Frequency (in Hz)'); 
ylabel('Data(f)');
grid on;
title('Single sided FFT plot for Original Data')

%% Filtering the data:
% perform band stop filter to remove the power interference.
f=600; %sampling frequency [Hz]% 
f_cutoff1=42.06;
f_cutoff2=42.1;
order=1;
fnorm1 = f_cutoff1/(f/2);
fnorm2 = f_cutoff2/(f/2);
[b1,a1] = butter(order,[fnorm1,fnorm2],'stop');
 y= filtfilt(b1,a1,p3);
 % perfrom the low pass filter to filter the high frequency noise.
 f_cutoff1=40;
f_cutoff2=0.2;
order=2;
fnorm1 = f_cutoff1/(f/2);
fnorm2 = f_cutoff2/(f/2);
[b1,a1] = butter(order,fnorm1,'low');
 y1= filtfilt(b1,a1,y);
 %% time domain analysis
 % select the range in which the peak detection window should perform
 select=y1(1:length(y1));
[ygh,x]=findpeaks(select,600,'Minpeakdistance',0.7);
% find the RR interval and HR
for i=1:1:size(x)-1
    xi(i)=x(i+1)-x(i);
    hr(i)=60/xi(i);
end
% perfrom the Statistical analysis for HRV
meanNNInterval=sum(xi)/(length(xi)*1000)
heartrate=60/(meanNNInterval*1000)
MinRange=min(xi)
MaxRange=max(xi)
DeviationinNN=std(xi)*1000
%% Frequency Domain Analysis
% to see the change in RR interval
k=1:1:length(xi)+1;

for k=1:1:length(xi)-1
    xii(k)=xi(k+1)-xi(k);
end

% resampling the RR interval data
rsf=20;
Xcs=resample(xi,rsf,1);

% see the change in RR interval
k=1:1:length(Xcs)+1;

for k=1:1:length(Xcs)-1
    xiii(k)=Xcs(k+1)-Xcs(k);
end
% plot the PSD
Fs=20;
% plot psd
N = length(Xcs); 
xdft = fft(Xcs);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) *abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(Xcs):Fs/2;


% PSD curve analysis in given frequency
yPower=10*log10(psdx);

ulf=abs(trapz(freq(2:10),yPower(2:10)))
lf=abs(trapz(freq(10:34),yPower(10:34)))
hf=abs(trapz(freq(34:91),yPower(34:91)))
tp=ulf+lf+hf
a=lf/hf
max=abs(min(yPower(10:60)))
totalpower=abs(trapz(freq(13:15),yPower(13:15)))
coherenceratio=totalpower/(tp-totalpower)

%% Plotting the graphs
% plot the original ECG data and Filtered ECG data vs Time
subplot(3,1,1)

plot(y11,p3)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Patient 1 Data Plot Vs Time');
grid on;
% subplot(3,2,2)

subplot(3,1,2)
plot(y11,y)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Filtered data: Bandstop Filter [42.06-42.1]Hz ');
grid on;
% subplot(3,2,4)
% plot(xfft,abs(fy1/max(fy1)));
subplot(3,1,3)
plot(y11,y1)
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Filtered data: low pass filter at 43Hz');
grid on;
figure()
% plot the Peak detetction of RR interval
findpeaks(select,600,'Minpeakdistance',0.7);
xlabel('Time (in Second)');
ylabel('Amplitude(v)');
title('Finding the peaks');
grid on
% plot the HR vs Time
figure()
plot(x(2:end),hr)
ylabel('Heart Rate');
xlabel(' Time (in seconds)');
title('HR vs Time');
figure()
% plot the rr interval vs time of resampled data and original data
subplot(2,1,1)
plot(x(2:end),xi)
xlabel('Time (in seconds)');
ylabel('R-R interval(x(i+1)-x(i))');
title('RR Interval plot');
subplot(2,1,2)
plot(1:1:length(Xcs),Xcs)
xlabel('Time (in seconds)');
ylabel('R-R interval Resampled plot at RSF=20Hz');
title('RR Interval resampled plot AT RSF=20Hz');
figure()
% plot the change in RR interval vs time for both original RR interval and
% change in RR interval
subplot(2,1,1)
plot(x(2:end-1),xii);   
xlabel('Time (in seconds)');
ylabel('R-R interval Change');
title('RR Interval Change plot');
subplot(2,1,2)
plot(1:1:length(Xcs)-1,xiii); 
xlabel('Time (in seconds)');
ylabel('R-R interval Change after Resampling & RSF=20Hz');
title('RR Interval Change  plot after resampling at 20Hz');
figure();
% plot the PSD curve.
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)');
end