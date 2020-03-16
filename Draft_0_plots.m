% This script is written on 5/13/19
% The purpose is to plot the figures use in the third article in my
% disseration. 
% A readme type discussion follows:

% ARforecastVals.mat: comes from aving the adjusted NMSE using the AR model in NMSE_t script
% ARXforecastVals.mat: comes from saving the adjusted NMSE using the AR model in NMSE_t script
% saveTime.mat: comes from saving the time from stakedAR_ARX.m
% FexConv.mat: comes from saving the true Fex from stackedARX_ARX.m with
% the correct settings
% FexEst: comes from saving the est Fex from stackedARX_ARX.m with
% the correct settings

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
plot(forecastVec,100.*ARnmse,'--',forecastVec,100.*ARXnmse)
legend('AR','ARX')
ylabel('Mean NMSE [%]')
xlabel('Prediction Horizon [s]')
title('Mean NMSE vs. Prediction Horizon')

c
%% Comparing ARX_AR vs AR only.
% Each one of the full NMSE matricies will be loaded. These are saved from
% the NMSE_t.m script 
load AR_NMSE.mat
AR_NMSE = NMSE;
load ARX_AR_NMSE
ARXAR_NMSE = NMSE;
Periods = [0,0,1.9,1.9,1.9,2.37,2.37,2.37,2.69,2.69,2.69,3.16,3.16,3.95,4.74,5.53,2.06,2.37];
deadTime = [0,0,2.53,2.53,2.53,2.03,2.03,2.03,1.79,1.79,1.79,1.52,1.52,1.22,1.01,0.89,2.33,2.03];
n = 20;
forecastVec = [1*Ts:Ts:n*Ts];
for k = 3:18 % create matrix of dimensionless forecast times.
dimForecast(:,k) = forecastVec/Periods(k);
end

extraHorizon = 0.5; % the time AR forecasts beyond dead time
cutTime = round((deadTime+extraHorizon)./Ts);
ARnans = AR_NMSE;
ARXARnans = ARXAR_NMSE;

for k = 1:length(forecastVec)
    for kk = 3:18
        if cutTime(kk) < forecastVec(k)/Ts
           ARnans(k:end,kk) = NaN;
           ARXARnans(k:end,kk) = NaN;
        end    
    end
    
end

for k = 1:length(forecastVec)
ARNMSE_adjusted(k) = nanmean(ARnans(k,3:18));
ARXARNMSE_adjusted(k) = nanmean(ARXARnans(k,3:18));
end


figure()
plot(forecastVec,100.*ARNMSE_adjusted,'b--',forecastVec,100.*ARXARNMSE_adjusted,'k')
xlabel('Prediction Horizon [s]')
ylabel('Mean NMSE [%]')
title('Mean NMSE, Horizon Extended 0.5 seconds')
legend('AR','ARX+AR')


%% Comparing AR and ARXAR
load ARNMSE_3_11.mat
AR9_NMSE = NMSE(:,3:11);
load ARXARNMSE_3_11.mat
ARXAR9_NMSE = NMSE(:,3:11);

Fs = 200/forecastParams.subNo;
Ts = 1/Fs;
n = length(AR9_NMSE);
forecastVec = [1*Ts:Ts:n*Ts];

for k = 2:3:10
    AR3_NMSE(:,round(k/3)) = mean(AR9_NMSE(:,k-1:k+1),2);
    ARXAR3_NMSE(:,round(k/3)) = mean(ARXAR9_NMSE(:,k-1:k+1),2);
end

figure()
Tp = 1.9;
subplot(3,1,1)
plot(forecastVec./Tp,100.*AR3_NMSE(:,1),'b--',forecastVec./Tp,100.*ARXAR3_NMSE(:,1),'r')
hold on
plot([1.33,1.33],[0,100],'k')
hold off
axis([0,forecastVec(end)./Tp,0,100])
ylabel('NMSE [%]')
title('Runs 3-5')

subplot(3,1,2)
Tp = 2.37;
plot(forecastVec./Tp,100.*AR3_NMSE(:,2),'b--',forecastVec./Tp,100.*ARXAR3_NMSE(:,2),'r')
hold on
plot([0.85,0.85],[0,100],'k')
hold off
axis([0,forecastVec(end)./Tp,0,100])
ylabel('NMSE [%]')
title('Runs 6-8')

subplot(3,1,3)
Tp = 2.69;
plot(forecastVec./Tp,100.*AR3_NMSE(:,3),'b--',forecastVec./Tp,100.*ARXAR3_NMSE(:,3),'r')
hold on
plot([0.66,0.66],[0,100],'k')
hold off
axis([0,forecastVec(end)./Tp,0,100])
ylabel('NMSE [%]')
title('Runs 9-11')
xlabel('Horizon/T_p')

% %% Testing individual runs
% close all
% for k = 1:9
%     figure(k)
%     plot(forecastVec,100.*AR9_NMSE(:,k),'b--',forecastVec,100.*ARXAR9_NMSE(:,k),'r')
% 
% end

%% 9 runs subplot to show individual wave climates.
load ARNMSE_10steps.mat
AR9 = NMSE(:,3:11);
load ARXNMSE_10steps.mat
ARX9 = NMSE(:,3:11);

Fs = 200/forecastParams.subNo;
Ts = 1/Fs;
n = 10;
forecastVec = [1*Ts:Ts:n*Ts];
deadTime = [0,0,2.53,2.53,2.53,2.03,2.03,2.03,1.79,1.79,1.79,1.52,1.52,1.22,1.01,0.89,2.33,2.03];
deadSteps = round(deadTime./Ts);
ARnans = 100.*AR9;
ARXnans = 100.*ARX9;

for k = 1:length(forecastVec)
    for kk = 1:9
        if deadSteps(kk+2) < forecastVec(k)/Ts % +2 adjusts for runNo index
           ARnans(k:end,kk) = NaN;
           ARXnans(k:end,kk) = NaN;
        end    
    end
    
end

figure()
subplot(3,1,1)
plot(forecastVec,mean(ARnans(:,1:3),2),'b--',forecastVec,mean(ARXnans(:,1:3),2),'r')
title('Mean NMSE for Runs 1-3, T_p=1.9s')
% axis([0,2.5,0,100])
ylabel('NMSE [%]')
legend('AR','ARX','Location','SouthWest')


subplot(3,1,2)
cutoffNum = 2;
plot(forecastVec(1:end-cutoffNum),mean(ARnans(1:end-cutoffNum,4:6),2),'b--',forecastVec(1:end-cutoffNum),mean(ARXnans(1:end-cutoffNum,4:6),2),'r')
title('Mean NMSE for Runs 4-6, T_p=2.37s')
% axis([0,2.5,0,100])
ylabel('NMSE [%]')
legend('AR','ARX','Location','SouthWest')


subplot(3,1,3)
cutoffNum = 3;
plot(forecastVec(1:end-cutoffNum),mean(ARnans(1:end-cutoffNum,7:9),2),'b--',forecastVec(1:end-cutoffNum),mean(ARXnans(1:end-cutoffNum,7:9),2),'r')
title('Mean NMSE for Runs 7-9, T_p=2.69s')
% axis([0,2.5,0,100])
xlabel('Prediction Horizon [s]')
ylabel('NMSE [%]')
legend('AR','ARX','Location','SouthWest')





