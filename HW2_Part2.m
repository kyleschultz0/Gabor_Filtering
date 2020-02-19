clc;clear all;close all

%% Reconstructing music score
%% Piano

tr_piano=16; % record time in seconds
y=audioread('music1.wav')'; Fs=length(y)/tr_piano;




t2=linspace(0,tr_piano,length(y)+1); t=t2(1:length(y)); 
k=(2*pi/tr_piano)*[0:length(y)/2-1 -length(y)/2:-1]; 
ks=fftshift(k);
tslide=0:0.1:tr_piano;
tau = 32;
max_freq=zeros(1,length(tslide));
amp = zeros(1,length(tslide));
Sgt_spec_piano = zeros(701440, length(tslide));



for j=1:length(tslide)
    g=exp(-tau(1)*(t-tslide(j)).^2); % Gabor
    Sg=g.*y;
    Sgt=fft(Sg);
    Sgt_spec_piano(:,j,1)=abs(fftshift(Sgt));
    [M,I]=max(abs(fftshift(Sgt)));
    max_freq(j)=-ks(I);
    amp(j) = M;
    %         figure(1)
    %         subplot(3,1,1), plot(t,y,'k',t,g,'r')
    %         subplot(3,1,2), plot(t,Sg,'k')
    %         subplot(3,1,3), plot(ks,Sgt_normal)
    %         xlim(10^3*[0 15])
    %         drawnow
end

amp_normal = amp/max(amp);
TF = islocalmin(amp_normal,'MinSeparation',0.2,'SamplePoints',tslide);
TF_shift = circshift(TF, 5);

E=329.63*ones(length(tslide));
D=293.66*ones(length(tslide));
C=261.63*ones(length(tslide));

figure(1);
plot(tslide, max_freq/(2*pi), tslide(TF_shift), max_freq(TF_shift)/(2*pi), 'r*');

legend("Music", "Key Strike")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
title("Time and Frequency of Notes Played")
xlim([0 15])
hold off

figure(2)
plot(tslide, amp_normal, tslide(TF), amp_normal(TF), 'r*')
title("Points of Minimum Amplitude") 
xlabel("Time (s)")
ylabel("Normalized Amplitude")
xlim([0 15])
% 
% for i=1:length(tau)
%     figure(3)
%     pcolor(tslide,ks,Sgt_spec_piano(:,:,i)), 
%     shading interp 
%     set(gca,'Ylim',10^3*[0 20],'Fontsize',[14]) 
%     title(sprintf('Piano, window width: %d',tau(i)))
%     colormap(hot)
% end

figure(3)
pcolor(tslide,ks/(2*pi),Sgt_spec_piano),
shading interp
set(gca,'Ylim',10^3*[0 1])
title("Piano Spectogram")
colormap(hot)
xlabel("Time (s)")
ylabel("Frequency (Hz)")

%% Recorder

tau = 32;
tr_rec=14; % record time in seconds
y=audioread('music2.wav')'; Fs=length(y)/tr_rec;


t2=linspace(0,tr_rec,length(y)+1); t=t2(1:length(y)); 
k=(2*pi/tr_rec)*[0:length(y)/2-1 -length(y)/2:-1]; 
ks=fftshift(k);
vt=fft(y);
tslide=0:0.02:tr_rec;
max_freq=zeros(1,length(tslide));
amp = zeros(1,length(tslide));
Sgt_spec_rec = zeros(length(t), length(tslide));



for i=1:length(tau)
    for j=1:length(tslide)
        g=exp(-tau(i)*(t-tslide(j)).^2); % Gabor
        Sg=g.*y;
        Sgt=fft(Sg);
        Sgt_spec_rec(:,j,i)=abs(fftshift(Sgt));
        Sgt_normal=abs(fftshift(Sgt));
        [M,I]=max(Sgt_normal);
        max_freq(j)=-ks(I);
        amp(j) = M;
%         figure(3)
%         subplot(3,1,1), plot(t,y,'k',t,g,'r')
%         subplot(3,1,2), plot(t,Sg,'k')
%         subplot(3,1,3), plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)))
%         drawnow
    end
end

amp_normal = amp/max(amp);
TF = islocalmin(amp_normal,'MinSeparation',0.2,'SamplePoints',tslide);
TF_shift = circshift(TF, 5);

figure(4)
plot(tslide, max_freq/(2*pi), tslide(TF_shift), max_freq(TF_shift)/(2*pi), 'r*');
legend('Signal', 'Blow in Recorder')
title("Time and Frequency of Notes Played")
xlabel("Time (s)")
ylabel("Max Frequency (Hz)")

figure(5)
plot(tslide, amp_normal, tslide(TF), amp_normal(TF), 'r*')
title("Points of Minimum Amplitude") 
xlabel("Time (s)")
ylabel("Normalized Amplitude")

% figure(6)
% pcolor(tslide,ks/(2*pi),Sgt_spec_rec),
% shading interp
% set(gca,'Ylim',10^3*[.5 1.5])
% title("Recorder Spectogram")
% colormap(hot)
% xlabel("Time (s)")
% ylabel("Frequency (Hz)")

% xlim([0 15])

% for i=1:length(tau)
%     figure(4)
%     pcolor(tslide,ks,Sgt_spec_rec(:,:,i)), 
%     shading interp 
%     set(gca,'Ylim',10^3*[0 20],'Fontsize',[14]) 
%     title(sprintf('Recorder, window width: %d',tau(i)))
%     colormap(hot)
% end