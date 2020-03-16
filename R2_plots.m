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

load AR_NMSE.mat
ARnans = NMSE;
load ARX_NMSE.mat
ARXnans = NMSE;

deadTime = [0,0,2.53,2.53,2.53,2.03,2.03,2.03,1.79,1.79,1.79,1.52,1.52,1.22,1.01,0.89,2.33,2.03];
deadSteps = round(deadTime./Ts);

for k = 1:length(forecastVec)
    for kk = 3:18
        if deadSteps(kk) < forecastVec(k)/Ts
            ARnans(k:end,kk) = NaN;
            ARXnans(k:end,kk) = NaN;
        end    
    end
    
end

for k = 1:length(forecastVec)
AR_std(k) = nanstd(ARnans(k,3:18));
ARX_std(k) = nanstd(ARXnans(k,3:18));
end


figure()
% plot(forecastVec,100.*ARnmse,'bo--',forecastVec,100.*ARXnmse,'ro-')
% boxplotData = [ARnans(1,3:end);ARnans(2,3:end);ARnans(3,3:end);ARnans(4,3:end);ARnans(5,3:end);ARnans(6,3:end);ARnans(7,3:end);ARnans(8,3:end);ARnans(9,3:end);ARnans(10,3:end);]'
% boxplot(boxplotData,forecastVec(1:10),'PlotStyle','compact')

hold on
plot(forecastVec,100.*ARnmse,'bo--','LineWidth',1)
plot(forecastVec,100.*ARXnmse,'ro-','LineWidth',1)
errorbar(forecastVec,100.*ARnmse,100.*AR_std,'--')
errorbar(forecastVec,100.*ARXnmse,100.*ARX_std)
% plot(forecastVec,100.*ARnmse+100.*AR_std,'c^-.','MarkerSize',4,'MarkerEdgeColor',[0,0,0.75])
% plot(forecastVec,100.*ARXnmse+100.*ARX_std,'r^')
% plot(forecastVec,100.*ARnmse-100.*AR_std,'cv-.')
% plot(forecastVec,100.*ARXnmse-100.*ARX_std,'rv')

hold off
legend('AR','ARX','location','SouthWest')
ylabel('Mean GOF [%]')
xlabel('Prediction Horizon [s]')
title('Mean GOF vs. Prediction Horizon')


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
AR_std(k) = nanstd(ARnans(k,3:18));
ARXAR_std(k) = nanstd(ARXARnans(k,3:18));
end

for k = 1:length(forecastVec)
ARNMSE_adjusted(k) = nanmean(ARnans(k,3:18));
ARXARNMSE_adjusted(k) = nanmean(ARXARnans(k,3:18));
end



figure()
hold on
plot(forecastVec,100.*ARNMSE_adjusted,'bo--','LineWidth',1)
plot(forecastVec,100.*ARXARNMSE_adjusted,'ko-','LineWidth',1)
errorbar(forecastVec,100.*ARNMSE_adjusted,100.*AR_std,'--')
errorbar(forecastVec,100.*ARXARNMSE_adjusted,100.*ARXAR_std,'k')
hold off
xlabel('Prediction Horizon [s]')
ylabel('Mean GOF [%]')
title('Mean GOF, Horizon Extended 0.5 seconds')
legend('AR','ARX+AR','Location','SouthWest')

% return
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
plot(forecastVec,mean(ARnans(:,1:3),2),'bo--',forecastVec,mean(ARXnans(:,1:3),2),'ro-')
% - std Dev
plot(forecastVec,mean(ARnans(:,1:3),2),'bo--',forecastVec,mean(ARXnans(:,1:3),2),'ro-')

title('Mean GOF for Runs 1-3, T_p=1.9s')
% axis([0,2.5,0,100])
ylabel('GOF [%]')
legend('AR','ARX','Location','SouthWest')


subplot(3,1,2)
cutoffNum = 2;
plot(forecastVec(1:end-cutoffNum),mean(ARnans(1:end-cutoffNum,4:6),2),'bo--',forecastVec(1:end-cutoffNum),mean(ARXnans(1:end-cutoffNum,4:6),2),'ro-')
title('Mean GOF for Runs 4-6, T_p=2.37s')
% axis([0,2.5,0,100])
ylabel('GOF [%]')
legend('AR','ARX','Location','SouthWest')


subplot(3,1,3)
cutoffNum = 3;
plot(forecastVec(1:end-cutoffNum),mean(ARnans(1:end-cutoffNum,7:9),2),'bo--',forecastVec(1:end-cutoffNum),mean(ARXnans(1:end-cutoffNum,7:9),2),'ro-')
title('Mean GOF for Runs 7-9, T_p=2.69s')
% axis([0,2.5,0,100])
xlabel('Prediction Horizon [s]')
ylabel('GOF [%]')
legend('AR','ARX','Location','SouthWest')





