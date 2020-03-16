% simKnownFexInput uses the convolution calculation of the excitation force
% to input a known input to the equations of motion (the function 'modeled response')
% to determine the process noise between the model and the actual
% experimental data.

clear all
close all
clc

% TempTimeDelays = [0.25:0.005:0.7];

rng(1)

% for jjj = 1:length(TempTimeDelays)

bodynumber = 2
for jj = bodynumber:bodynumber % for each of the bodies

for j = 7:7
% for j = 7:7
%% Initialization and Data Generation
% generate path to needed folders
Paths.exampleDir = cd;
Paths.functionsDir = fullfile(Paths.exampleDir, '..', filesep,'..', filesep, 'Functions'); 
Paths.dataDir = fullfile(Paths.exampleDir, '..', filesep,'..', filesep, 'Data'); 
% Paths.functionsDir = fullfile(Paths.exampleDir, '..', filesep, 'Functions'); % original path for Dev/Test
Paths.outputDir = fullfile(Paths.exampleDir, filesep, 'Outputs'); % original path for Dev/Test
% Paths.dataDir = fullfile(Paths.exampleDir, '..', filesep, 'Data'); % original path for Dev/Test
% generate GLOBAL VARIABLES
global Param Time Func Data SpectrumData bodyData
Param.bodyNumber = jj;


% Simulation and data choices
name.body = {'FPS','OWC'};
model = name.body(jj);
bodyFunc = @generateBody_BosmaFPS; % this should be unique to the application
% bodyString: decide whether you want to run NEMOH
bodyString = 'check'; % options: skip, run, or check
% dataString decide wheather you want to load measurements bosmaData or
% simulate
dataString = 'check'; %  options: read, check, or sim
% measurements are the values evaluated by the EKF. These can be
% experiemental or numerical
savedDataName = ' '; %'RefSim_JP1.mat';
measString = 'loadExper'; % options: loadSim, loadExper, sim
% decide if you want to generate the functions
genFunString = 'gen'; % options: gen, check
%

waveBool = [1]; % options: 1 -> generate wave, 0-> no wave, 2 -> regular wave with unit amplitude
mooringBool = [0]; % 0-> taught & 1-> catenary
RTSBool = [0]; % 0-> Do NOT implement RTS, 1-> run RTS Smoother
% this is used in the estimateVel function
Param.velEstimateVal = 1; % 0-> uses particleDynamics_spec 1-> uses method 1 2-> method 2...

% generate MATLAB Path
% addedPath = genpath('..\Functions');
addpath(Paths.functionsDir)
addpath('Outputs')
[status, msg, msgID] = mkdir('Outputs');


%% Enter Parameters
Param.numBodies = length(name.body);
Param.RTSBool = RTSBool;

% Choosing specific data
% runNo = '125'; Param.wavePeriod = 2.69; Param.waveHs = 0.3;
% Choosing specific data
% Param.runNo = '125'; Param.wavePeriod = 2.69; Param.waveHs = 0.3;
% Param.runNo = '113'; Param.wavePeriod = 3; Param.waveHs = 0.2; % regular wave
% Param.runNo = '099'; Param.wavePeriod = 2.06; Param.waveHs = 0.1; % regular wave
% Param.runNo = '117'; Param.wavePeriod = 4.74; Param.waveHs = 0.2; % regular wave
% Param.runNo = '162'; Param.waveperiod = 3; Param.waveHs = 0.2; % regular wave O ring
% Param.runNo = '176'; Param.wavePeriod = 2.69; Param.waveHs = 0.3; % regular wave with O-ring
% Param.runNo = '122'; Param.wavePeriod = 2.21; Param.waveHs = 0.3; % regular TFI 30%
% Param.runNo = '173'; Param.wavePeriod = 2.21; Param.waveHs = 0.3; % regular O-ring
% 
% oRingIndex = 18;
oRingIndex = j
waveRuns = [1,180,0.075,1.1;
            2,181,0.175,1.42;
            3,182,0.05,1.9;
            4,183,0.175,1.9;
            5,184,0.3,1.9;
            6,185,0.050,2.37;
            7,186,0.175,2.37;
            8,199,0.3,2.37;
            9,189,0.05,2.69;
            10,190,0.175,2.69;
            11,191,0.3,2.69;
            12,192,0.05,3.16;
            13,193,0.175,3.16;
            14,194,0.175,3.95;
            15,195,0.175,4.74;
            16,196,0.175,5.53;
            17,197,0.375,2.06;
            18,198,0.375,2.37];
Param.runNo = num2str(waveRuns(oRingIndex,2)); Param.wavePeriod = waveRuns(oRingIndex,4); Param.waveHs = waveRuns(oRingIndex,3);

% This is the result of the wave sensor and the OWC being several feet
% apart.
Param.oRingIndex = oRingIndex;
% Param.timeDelay = ones(1,18)*TempTimeDelays(jjj);
% Param.timeDelay = [0.65,0.45,0.55,0.4,0.30,0.5,0.4,0.35,0.5,0.40,0.35,0.45,0.4,0.4,0.45,0.45,0.3,0.3];
Param.timeDelay = [0.6550,0.4400,0.530,0.4150,0.280,0.5050,0.4200,0.3300,0.4950,0.420,0.3350,0.4750,0.42,0.42,0.44,0.44,0.2850,0.3050];

name.sensorFile = strcat('Run_', Param.runNo, '.csv');
name.waveFile = strcat('Run_', Param.runNo, '.txt');
name.motionFile = strcat('Run_', Param.runNo, '_6D.tsv'); 
name.runNo = Param.runNo;


Param.samplingFreq.motion = 200;
Param.samplingFreq.wave = 128;
Param.samplingFreq.loadCell = 128;
Param.peakFrequency = 1/Param.wavePeriod;
Param.waveLength = 9.81*Param.wavePeriod/(2*pi);

Time.subSampleNo = 1; % uses one out of every subSampleNo
Time.dT = Time.subSampleNo/Param.samplingFreq.motion;
Time.simStepsPerDT = 50;
Time.dtSim = Time.dT/Time.simStepsPerDT;
Time.start = 245; % min 1 period into timespan
Time.end = 260;% max 1 period from the end of the timespan
Time.t = [Time.start:Time.dT:Time.end];

% Number of states and control inputs
Param.sizeRadSS = 4;
Param.numIDStates = 1;
Param.numStates = 2+Param.numIDStates+Param.sizeRadSS; % 2 for dynamics 
Param.numOutputs = 2; 
Param.numNonRad = Param.numStates-Param.sizeRadSS;

% initialize EstStates
% Both of these global variables are modified during the smoother_RTS funciton
Est.states = zeros(Param.numStates,length(Time.t));
Est.vel = zeros(1,length(Time.t)); % updated in the estimateVel function
SmoothCounter = 0;

% Param structure data Properties of the experiment
Param.rho = 1000;
Param.g = 9.81;
Param.mass(1) = 9.582; % 1/10 scale FPS
Param.mass(2) = 76.294; % 1/10 scale OWC
Param.upperPanelLength(1) = 0.51/(sqrt(2)+1); % originally this was 0.50
Param.upperPanelLength(2) = 0.75/(sqrt(2)+1);
Param.lowerPanelLength(1) = 0.35/(sqrt(2)+1);
Param.lowerPanelLength(2) = 0.50/(sqrt(2)+1);
Param.avgPanelLength(1) = (Param.upperPanelLength(1)+Param.lowerPanelLength(1))/2;
Param.avgPanelLength(2) = (Param.upperPanelLength(2)+Param.lowerPanelLength(2))/2;
Param.lowerArea(1) = 2*(1+sqrt(2))*Param.lowerPanelLength(1)^2;
Param.upperArea(1) = 2*(1+sqrt(2))*Param.upperPanelLength(1)^2;
Param.avgArea(1) = 2*(1+sqrt(2))*Param.avgPanelLength(1)^2;
Param.lowerArea(2) = 2*(1+sqrt(2))*Param.lowerPanelLength(2)^2;
Param.upperArea(2) = 2*(1+sqrt(2))*Param.upperPanelLength(2)^2-Param.lowerArea(2);

Param.buoyDepth = 0.075; % specified as depth from surface
Param.waveNo = Param.peakFrequency*2*pi; % deep water assumption implicit
Param.attenFactor = exp(-Param.buoyDepth*Param.waveNo); % calculated on peak frequency
Param.attenFactor = 1; % this removes the attenfactor influence while keeping it for further implementation

Param.addedMass(1) = 21.15;
% Param.addedMass(2) = 29;
Param.addedMass(2) = 0;
% Param.exciteTerm = 1800;
Param.weight(1) = Param.mass(1)*Param.g;
Param.weight(2) = Param.mass(2)*Param.g;
Param.totalInertia(1) = Param.mass(1)+Param.addedMass(1);
Param.totalInertia(2) = Param.mass(2)+Param.addedMass(2);
% Param.hydrostaticCoeff = Param.avgArea*Param.rho*Param.g;
Param.hydrostaticCoeff(1) = Param.upperArea(1)*Param.rho*Param.g;
Param.hydrostaticCoeff(2) = Param.upperArea(2)*Param.rho*Param.g;
Param.heaveMoorStiff(1) = 200; % same moorings used for both buoys
Param.heaveMoorStiff(2) = 150;

% Param.dispVol = 11322077e-9; % from solidworks for the FPS
% Param.Delta = Param.dispVol*Param.rho;
% Param.dragFactor(1) = 1/2*Param.rho*Param.upperArea(1);
Param.C_D = 0.4;
Param.C_D(2) = Param.C_D(1);
Param.dragFactor(1) = Param.C_D(1)*1/2*Param.rho*Param.upperArea(1);
Param.dragFactor(2) = Param.C_D(2)*1/2*Param.rho*Param.upperArea(2);
Param.linearDrag(1) = 60;
Param.linearDrag(2) = 120;


% Initial Condtions 
Param.IC.est = [0.1, 0.1, 1,0,0,0,0];
Param.IC.true = [0, 0, 0.4, 0,0,0,0]; % parameters used when simulating 'data'
% Param.IC.P = 0.1*eye(Param.numStates);
Param.IC.P = zeros(Param.numStates);
Param.IC.P(1,1) = 0.1; Param.IC.P(2,2) = 0.1; Param.IC.P(3,3) = 0.1;
% Param.IC.P = 1*eye(Param.numStates)

% Artificial noise addition
Param.Ra = [0,0]; % artificial noise added as measurement noise
Param.Ga = [0,0]; % artificial noise added as 'process noise'


%% Data initializations ___________________________________________________
% Generating or loading body data
[bodyData] = check_generateBody(bodyString,name,bodyFunc,Paths,Param.sizeRadSS);
for k = 1:length(name.body)
    Param.body(k).A_r = bodyData.body(k).A;
%     Param.body(k).A_r = zeros(4,4);
    Param.body(k).B_r = bodyData.body(k).B;
%     Param.body(k).B_r = zeros(4,1);
    Param.body(k).C_r = bodyData.body(k).C;
%     Param.body(k).C_r = zeros(1,4);
%     Param.body(k).D_r = bodyData.body(k).D;
    Param.body(k).D_r = 0; % see Perez and Fossen body of work
end
% Generating or loading Bosma Data
[Data] = check_readBosmaData(dataString,name,Paths);
% Defining simulated spectrum data
[SpectrumData] = generateSpectrum(waveBool,bodyData);
% ___________________________________________________________________

% EKF initializations
% Param.options = odeset();
% Param.R = (0.0006075)^2*eye(Param.numOutputs); % variance of the residuals
% display('modified R')
% process noise covariance matrix
% Param.GQGt = zeros(Param.numStates);

% determining the input wave variance
%     for k = 1:length(Time.t)
%         [pos,vel,acc] = particleDynamics_spec(Time.t(k));
%         incidentWaveTemp(k) = pos;
%     end
% stdDevInputPos = std(incidentWaveTemp);

[savedData] = check_measurements(measString,savedDataName,Paths,[Time.start,Time.end]);
Measurements = savedData;

%%

% Generating simulation based on convolution calculated force
% T = [Time.start:1/Time.simStepsPerDT:Time.end];
q0 = zeros(1,6);

% if jj ==1
[t,q] = ode45('modeledResponse2bodies',Time.t,q0); 

% selecting a subset of the data
% Measurements = Measurements;
% q = q(1:Time.subSampleNo:end,:);

% queing up a bunch of data
if Param.bodyNumber == 1
for k = 1:length(Time.t)
    U1(k,1) = interp1(Data.wave.time,Data.wave.incidentFPS,Time.t(k));
    U2(k,1) = interp1(Data.wave.time(1:end-1),Data.wave.incidentFPSVel,Time.t(k));
    U3(k,1) = interp1(Data.wave.time(1:end-2),Data.wave.incidentFPSAcc,Time.t(k));    
end  
elseif Param.bodyNumber == 2
    for k = 1:length(Time.t)
    U1(k,1) = interp1(Data.wave.time,Data.wave.incidentOWC,Time.t(k));
    U2(k,1) = interp1(Data.wave.time(1:end-1),Data.wave.incidentOWCVel,Time.t(k));
    U3(k,1) = interp1(Data.wave.time(1:end-2),Data.wave.incidentOWCAcc,Time.t(k));    
    end
end
% coefficient of determination for position
num = sum((Measurements(:,1)-q(:,1)).^2);
den = sum((Measurements(:,1)).^2);
NMSE(j) = 1-(num/den)


%% making a data structure to load wave information
runInfo(j).Hs = Param.waveHs;
runInfo(j).Hs_ScaledUp = 10*Param.waveHs;
runInfo(j).Ts = Param.wavePeriod;
runInfo(j).Ts_ScaledUp = sqrt(10)*Param.wavePeriod;
runInfo(j).waveLength = (2*pi)/Param.waveNo;
runInfo(j).steep = Param.waveHs/runInfo(j).waveLength;
runInfo(j).sizeRatio = 0.51/runInfo(j).waveLength;

posError = Measurements(:,1)-q(:,1);
velError = Measurements(:,2)-q(:,2);

if Param.bodyNumber == 1
Disp = Data.wave.incidentFPS;
DispTime = Data.wave.time;
F_exIRF = bodyData.body(Param.bodyNumber).FIRF;
IRFtime = bodyData.body(Param.bodyNumber).irfTime;
elseif Param.bodyNumber == 2
Disp = Data.wave.incidentOWC;
DispTime = Data.wave.time-Param.timeDelay(Param.oRingIndex);
F_exIRF = bodyData.body(Param.bodyNumber).FIRF;
IRFtime = bodyData.body(Param.bodyNumber).irfTime;
end

for k = 1:length(Time.t)

    exforce(k) = excitationConv(Time.t(k),F_exIRF,IRFtime,Disp,DispTime);
%     exforce_sim(k) = excitationForce_spec(Time.t(k));
end

% runInfo(j).FexSTD = std(exforce_sim); % carefull attention to the df vs dW needs to be paid here
runInfo(j).FexSTD = std(exforce);
runInfo(j).posNoiseSTD = std(posError);
runInfo(j).velNoiseSTD = std(velError);
runInfo(j).NMSE = NMSE(j);

% TempNMSE(j,jjj) = NMSE(j);
% saving convolution calculation of excitation force
Excitation.Fex = exforce;

save(string(strcat('Fex_Conv_',model,'_Run',num2str(j),'_Time',num2str(Time.start),'_',num2str(Time.end),'.mat')),'Excitation')


end
%% plotting results
figure()
plot(Time.t,q(:,1),Time.t,Measurements(:,1))
legend('Sim','Exper')
title('Position')
xlabel('s')
ylabel('m')

figure()
plot(Time.t,q(:,2),Time.t,Measurements(:,2))
legend('Sim','Exper')
title('Velocity')
xlabel('s')
ylabel('m/s')

% error in position
figure()
error = Measurements(:,1)-q(:,1);
meanMeas = mean(Measurements(:,1))
meanSim = mean(q(:,1))
plot(Time.t,error)
title('error')
xlabel('s')
ylabel('m')

% error in velocity
figure()
errorVel = Measurements(:,2)-q(:,2);
plot(Time.t,errorVel)
title('error')
xlabel('s')
ylabel('m/s')


% autocorr(error)
PosVar = var(error)
PosSTD = std(error)

VelVar = var(errorVel)
VelSTD = std(errorVel)

figure()
plot(Time.t,q(:,1).*2000,Time.t,2000.*Measurements(:,1),Time.t,exforce)
legend('simPos*2000','measPos*2000','MeasFex')

% figure()
% plot(Data.wave.time,Data.wave.incidentFPS,Time.t,exforce./1000)
% legend('wave','excitation')
% axis([Time.start,Time.end,-1,1])

save(string(strcat('runStruct',model,'.mat')),'runInfo')
% save('runStructFPS.mat','runInfo')
% save('runStructSub2.mat','runInfo')

end

% end

