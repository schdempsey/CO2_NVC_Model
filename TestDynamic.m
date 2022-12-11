addpath(strcat(pwd,'\Cost_Functions'));
addpath(strcat(pwd,'\Exp_Data'));
addpath(strcat(pwd,'\Model_Files'));
addpath(strcat(pwd,'\Other_Scripts'));
Constants=load('Constants.mat');
Con=Constants.Constants;
load('theta_Shmuel.mat')
kCO2=log10(0.02);
Dng=log10(50);
Dgc=log10(50);
kCx26g=log10(0.001);
sinkATP=log10(1);
kATP=log10(100000);
kCx261=log10(0.5);
kCx262=log10(0.5);
theta=[X(1:37)';
    kCO2;
    Dng;
    Dgc;
    kCx26g;
    sinkATP;
    kATP;
    kCx261;
    kCx262];
options = amioption('sensi',0,...
    'maxsteps',1e4);
options.sensi = 0;
options.nmaxevent = 1e3;
% Steady state simulation
Ca_start = 10;
cCO2in=0.3; %percent CO2 in
sol = simulate_SSM_H22_A3(inf,theta(4:43),[Ca_start,Con,cCO2in],[],options);
HbO_0 = sol.y(2);
HbR_0 = sol.y(3);
SaO2_0 = sol.y(4);
ScO2_0 = sol.y(5);
SvO2_0 = sol.y(6);
options.x0 = sol.x(end,:).';
TE = 20*10^-3;       B0 = 7;
stimend=[0 20 40 14 10];
Constants = [sol.x(end,[11 9 13 24]), Ca_start, stimend(3)-stimend(2),...
    Con, HbO_0, HbR_0, SaO2_0, ScO2_0, SvO2_0, TE, B0, cCO2in];
% alter simulation tolerances, DAE solver can not handle the default values
options.atol = 1e-6;
options.rtol = 1e-12;
solReal1 = simulate_Dyn_H22_A3(0:1:50,[theta(4:45) ;theta(1:3)],Constants,[],options);

figure(1)
subplot(1,2,2)
plot(solReal1.t+30,solReal1.y(:,6),'r')
%R1=plot(solReal2.t+30,solReal2.y(:,6),'r');
%R1.Color = [R1.Color 0.4];
subplot(1,2,2)
plot(solReal1.t+30,solReal1.y(:,5),'r')    
%R1=plot(solReal2.t+30,solReal2.y(:,5),'r');
%R1.Color = [R1.Color 0.4];






