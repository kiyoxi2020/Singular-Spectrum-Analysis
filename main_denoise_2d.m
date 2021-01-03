
clear
close all
%% Load data
data = load('data.mat');
data_noisy = data.data_noisy;
data_clean = data.data_clean;
%% denoising
flow=0;fhigh=250;dt=0.004;N=3;
data_denoise=fxyssa_denoise_2d(data_noisy,flow,fhigh,dt,N);

figure;
imagesc([data_clean(:,:,1),...
    data_noisy(:,:,1),...
    data_denoise(:,:,1),...
    data_denoise(:,:,1)-data_noisy(:,:,1)]);
caxis([-0.5,0.5]);colormap gray

