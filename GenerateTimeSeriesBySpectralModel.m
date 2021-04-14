% -----------------------------
% Script: Generate Time Series By Spectral Model
% ------------
% History:
% v01:  By Feng Guo on 17-Mar-2021 supported by David Schlipf
% Flensburg University of Applied Sciences
% This research has received funding from the European Union's Horizon 2020 research and innovation program 
% under the Marie Sk{\l}odowska-Curie grant agreement No. 858358 (LIKE -- Lidar Knowledge Europe).
% ----------------------------------
clear all
close all
clc

load('ParaFit.mat')


%%% simulate by frequency and rocof spectra------------------------
%% Configuration 
dt          = 0.0625;         % time step
f_max       = 1/dt/2;         % Nyquist frequency
T           = 600;            % Time length
df          = 1/T;            % frequency step
f           = df:df:f_max;    % frequency vector
n_f         = length(f);      % number of frequency bins
t           = 0:dt:T-dt;      % time vector
n_t         = length(t);      % number of time steps


%% illustration of fitted function
a1          = x_opt(1);
b1          = x_opt(2);
a2          = x_opt(3);
b2          = x_opt(4);
f1          = x_opt(5);
f2          = x_opt(6);

fc          = sqrt(f1*f2);
r           = f2/f1;
t_b         = @(f) exp(-(f/fc).^(2/log10(r)));
S1          = @(f) a1*f.^b1;
S2          = @(f) a2*f.^b2;
S_b         = @(f) S1(f).* t_b(f) + S2(f).* (1-t_b(f));
S_f         = S_b(f);
S_rocof     =  (4*pi^2.*f.^2.*S_f)';
S_f         = S_f';


figure(1)

subplot(4,1,1)
grid on; hold on; box on;
plot(f,S_f)
set(gca,'xscale','log')
xlabel('frequency [Hz]')
ylabel('PSD [Hz^2/Hz]')
title('frequency PSD')

subplot(4,1,2)
plot(f,S_rocof)
grid on; hold on; box on;
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('frequency [Hz]')
ylabel('PSD [(Hz/s)^2/Hz]')
title('RoCoF PSD')
%% simulate time series
% init RNG
rng(1);

std_f_model   = sqrt(sum(S_f)*df);
std_f_rocof   = sqrt(sum(S_rocof)*df);

scale_f       = 0.2;      % should be adjusted based on measured data 
scale_rocof   = 0.05;     % should be adjusted based on measured data

% generation
n_FFT         = n_t;
A_rocof       = sqrt(2*S_rocof*df)*scale_rocof;
A_f           = sqrt(2*S_f*df)*scale_f;

% phase angles
Phi           = rand(1,n_f)*2*pi;
Phi2          = rand(1,n_f)*2*pi;

fgrid         = n_f*[0 A_f'].*exp([0,1i*Phi]);
Rocof         = n_f*[0 A_rocof'].*exp([0,1i*Phi]);

f_var         = ifft(fgrid,n_t,'symmetric');
rocof         = ifft(Rocof,n_t,'symmetric');


subplot(4,1,3)
grid on; hold on; box on;
plot(t,f_var)
xlabel('Tims [s]')
ylabel('Hz')
title('frequency time series')

subplot(4,1,4)
grid on; hold on; box on;
plot(t,rocof)
xlabel('Tims [s]')
ylabel('Hz/s')
title('RoCoF time series')
