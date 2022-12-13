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
theta=[X(1:37)';
    kCO2;
    Dng;
    Dgc;
    kCx26g;
    sinkATP;
    kATP];
options = amioption('sensi',0,...
    'maxsteps',1e4);
options.sensi = 0;
options.nmaxevent = 1e3;
% Steady state simulation
Ca_start = 10;
cCO2in=0.3;
% steady state simulation
sol = simulate_SSM_H22_A3(inf,theta(4:end),[Ca_start,Con,cCO2in],[],options);





