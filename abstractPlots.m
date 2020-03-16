% This function is written on 4/26/19 and is intended to take data saved
% from the stackedAR_ARX.m function. 

clear all
close all
clc
format compact
addpath(fullfile(cd, '..', filesep, 'Functions'))

load ARforecastVals.mat
    ARforecast = saveVals;
load ARXforecastVals.mat
    ARXforecast = saveVals;
load saveTime.mat
    Time = saveTime;
load FexConv.mat 
    Conv = saveConv;
load FexEst.mat
    Est = saveFex;

figure()
hold on
plot(Time,Est,'k')
plot(Time,ARforecast,'b--')
plot(Time,ARXforecast,'r-.')
legend('True','AR','ARX','Location','NorthWest')
title('2 Second Prediction Comparison')
ylabel('Excitation Force [N]')
xlabel('Time [s]')
axis([352,392,-330,330])


%% NMSE(t) plot
forecastParams.subNo = 50;
Fs = 200/forecastParams.subNo;
Ts = 1/Fs;
n = 10;
forecastVec = [1*Ts:Ts:n*Ts];

load ARNMSE_adjusted.mat
ARnmse = NMSE_adjusted;
load ARXNMSE_adjusted.mat
ARXnmse = NMSE_adjusted;

figure()
plot(forecastVec,ARnmse,'--',forecastVec,ARXnmse)
legend('AR','ARX')
ylabel('NMSE [%]')
xlabel('Prediction Horizon [s]')
title('NMSE vs. Prediction Horizon')














