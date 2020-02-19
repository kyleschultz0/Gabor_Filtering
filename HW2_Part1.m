clc; clear all; close all

%% Part 1

load handel
v = y'/2;
L=length(v)/Fs;
t2=linspace(0,L,length(v)+1); t=t2(1:length(v)); 
k=(2*pi/L)*[0:length(v)/2 -length(v)/2:-1]; 
ks=fftshift(k);
vt=fft(v);

p8 = audioplayer(v,Fs);
playblocking(p8)

tslide=0:0.1:L;
% tau = [1 8 64];
tau=64;
Sgt_spec_gauss=[];

Sp=spectrogram(v);
figure(1)
pcolor(abs(Sp))
shading interp
set(gca,'Fontsize',[14]) 
% set(gca,'Ylim',10^3*[0 20],'Fontsize',[14]) 
colormap(hot)

%% Effect of varying window width
for i=1:length(tau)
    for j=1:length(tslide)
        g=exp(-tau(i)*(t-tslide(j)).^2); % Gabor 
        Sg=g.*v; 
        Sgt=fft(Sg); 
        Sgt_spec_gauss(:,j,i)=abs(fftshift(Sgt)); 
        subplot(3,1,1), plot(t,v,'k',t,g,'r'), title("Sliding Filter Over Signal"), xlabel("Time")
        subplot(3,1,2), plot(t,Sg,'k'), title("Windowed Signal"), xlabel("Time")
 	    subplot(3,1,3), plot(ks,abs(fftshift(Sgt))/max(abs(Sgt))), title("Windowed Spectrum"), xlabel("Frequency")
        drawnow
        pause(1)
    end
end

figure(2)
for i=1:length(tau)
%     subplot(3,1,i)
    pcolor(tslide,ks,Sgt_spec_gauss(:,:,i)), 
    shading interp 
    set(gca,'Ylim',10^3*[0 5]) 
    title("Spectogram")
    colormap(hot)
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
end

%% Varying sampling frequency

t_sample = [0.01 0.1 2];
tau = 8;
Sgt_spec=[];
for i = 1:length(t_sample)
    tslide = 0:t_sample(i):L;
        for j=1:length(tslide)
            g=exp(-tau*(t-tslide(j)).^2); % Gabor 
            Sg=g.*v; 
            Sgt=fft(Sg); 
            Sgt_spec(:,j,i)=abs(fftshift(Sgt)); 
        end
end

t_sample_label = [10 100 2000];

figure(3)
sgtitle("Varying Sample Time", 'FontSize' ,14,"FontWeight","bold")
for i=1:length(t_sample)
    tslide = 0:t_sample(i):L;
    subplot(3,1,i)
    pcolor(tslide,ks,Sgt_spec(:,1:length(tslide),i)), 
    shading interp 
    set(gca,'Ylim',10^3*[0 5],'Fontsize',10) 
    title(sprintf('Sampling time: %d ms',t_sample_label(i)), "FontSize",12,"FontWeight","normal")
    colormap(hot)
    xlabel("Time (s)",'Fontsize',11)
    ylabel("Frequency (Hz)",'Fontsize',11)
end

tslide=0:0.1:L;

%% Mexican hat wavelet

% sigma = [4 1 0.5 0.25 0.1 0.05];
sigma=0.0565;
for i=1:length(sigma)
    for j=1:length(tslide)
        g_sobrero=2/(sqrt(3*sigma(i)).*pi.^(1/4)).*(1-((t-tslide(j))/sigma(i)).^2)...
            .*exp(-(t-tslide(j)).^2/(2.*sigma(i).^2)); % Mexican hat wavelet 
        Sg=g_sobrero.*v; 
        Sgt=fft(Sg); 
        Sgt_spec_hat(:,j,i)=abs(fftshift(Sgt)); 
%         figure(4)
%         subplot(3,1,1), plot(t,v,'k',t,g_sobrero,'r')
%         subplot(3,1,2), plot(t,Sg,'k')
%  	    subplot(3,1,3), plot(ks,abs(fftshift(Sgt))/max(abs(Sgt))) 
%         drawnow
    end
end

% for i=1:length(sigma)
%     figure(5)
%     subplot(3,3,i)
%     pcolor(tslide,ks,Sgt_spec_hat(:,:,i)), 
%     shading interp 
%     set(gca,'Ylim',10^3*[0 20],'Fontsize',[14]) 
%     title(sprintf('Window S.D.: %d s',sigma(i)))
%     colormap(hot)
% end

%% Shannon (step) filter

% t_width=[0.125 0.5 1];
t_width = 0.5;
width=t_width.*Fs;
for i=1:length(t_width)
    n=length(tslide);
    tslide = tslide+t_width(i)/2;
    for j=1:n
        g_shannon=zeros(1,length(v));
        g_shannon((tslide(j)*Fs-width(i)/2+1):1:(tslide(j)*Fs+width(i)/2+1))...
          =ones(1,width(i)+1);
        g_shannon=g_shannon(1:length(v));
        Sg=g_shannon.*v; 
        Sgt=fft(Sg); 
        Sgt_spec_shannon(:,j,i)=abs(fftshift(Sgt)); 
%         figure(6)
%         subplot(3,1,1), plot(t,v,'k',t,g_shannon,'r')
%         subplot(3,1,2), plot(t,Sg,'k')
%         subplot(3,1,3), plot(ks,abs(fftshift(Sgt))/max(abs(Sgt))) 
%         drawnow
    end
end

% for i=1:length(t_width)
%     figure(7)
%     subplot(3,1,i)
%     pcolor(tslide,ks,Sgt_spec_shannon(:,:,i)), 
%     shading interp 
%     set(gca,'Ylim',10^3*[0 20],'Fontsize',[14]) 
%     title(sprintf('Window Width: %d s',t_width(i)))
%     colormap(hot)
% end

%% COMPARING SPECTROGRAM OF DIFFERENT FILTERS

figure(8)
sgtitle("Comparing Different Filters", "FontSize",14,"FontWeight","bold")
subplot(3,1,1)
pcolor(tslide,ks,Sgt_spec_gauss(:,:,1)),
shading interp
set(gca,'Ylim',10^3*[0 5],'Fontsize',[10])
title("Gaussian", "FontSize",12,"FontWeight","normal")
xlabel("Time (s)",'Fontsize',11)
ylabel("Frequency (Hz)",'Fontsize',11)
subplot(3,1,2)
pcolor(tslide,ks,Sgt_spec_hat(:,:,1)),
shading interp
set(gca,'Ylim',10^3*[0 5],'Fontsize',[10])
title("Mexican Hat", "FontSize",12,"FontWeight","normal")
xlabel("Time (s)",'Fontsize',11)
ylabel("Frequency (Hz)",'Fontsize',11)
subplot(3,1,3)
pcolor(tslide,ks,Sgt_spec_shannon(:,:,1)),
shading interp
set(gca,'Ylim',10^3*[0 5],'Fontsize',[10])
title("Shannon", "FontSize",12,"FontWeight","normal")
xlabel("Time (s)",'Fontsize',11)
ylabel("Frequency (Hz)",'Fontsize',11)
colormap(hot)

%% DEV ONLY

tau=12.56;
sigma=0.0565;
t_width=0.5;
width=t_width.*Fs;
tslide=4.5;

g=exp(-tau*(t-tslide).^2); % Gabor

g_sobrero=2/(sqrt(3*sigma).*pi.^(1/4)).*(1-((t-tslide)/sigma).^2)...
    .*exp(-(t-tslide).^2/(2.*sigma.^2)); % Mexican Hat

g_shannon=zeros(1,length(v));
g_shannon((tslide*Fs-width/2+1):1:(tslide*Fs+width/2+1))...
    =ones(1,width+1);
g_shannon=g_shannon(1:length(v));

figure(9)
plot(t,g,t,g_sobrero,t,g_shannon)
title("Different Filters in Time Domain")
legend("Gaussian", "Mexican Hat", "Shannon")
xlabel("Time")
ylabel("Amplitude")
xlim([2 7])




