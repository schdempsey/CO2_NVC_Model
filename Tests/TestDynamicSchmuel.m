clear all;clc;
addpath(strcat(pwd,'\Cost_Functions'));
addpath(strcat(pwd,'\Exp_Data'));
addpath(strcat(pwd,'\Model_Files'));
addpath(strcat(pwd,'\Other_Scripts'));
Constants=load('Constants.mat');
Con=Constants.Constants;
load('theta_Shmuel.mat')
kCO2   =log10(0.02);
Dng    =log10(50);
Dgc    =log10(50);
kCx26g =log10(0.001);
sinkATP=log10(1);
kATP   =log10(100000);
kCx261 =log10(0.5);
kCx262 =log10(0.5);
newParams=[kCO2;
    Dng;
    Dgc;
    kCx26g;
    sinkATP;
    kATP;
    kCx261;
    kCx262];
theta=[X(1:37)';newParams];
options = amioption('sensi',0,...
    'maxsteps',1e4);
options.sensi = 0;
options.nmaxevent = 1e3;
% Steady state simulation
Ca_start = 10;
cCO2in=0.3; %percent CO2 in

sol = simulate_SSM_H22_A3(inf,theta(4:end),[Ca_start,Con,cCO2in],[],options);
HbO_0 = sol.y(2);
HbR_0 = sol.y(3);
SaO2_0 = sol.y(4);
ScO2_0 = sol.y(5);
SvO2_0 = sol.y(6);
%cCO2c0 = sol.y(7);
%%
options.x0 = [0;0;0; %Activities
    0.00296945972646950; %4 Ca NO
    0.00671011508027784; %Ca NPY
    0.114498553469019; %Calcium Pyr Changed
    1.75110047557811e-09;% 7AA
    0.531373632118709;% 8 PGE
    0.000581299567725387;% 9 PGE vsm
    0.0127733784689748;%10 NO
    7.35304240575339e-06;%11 NO vsm
    0.00175808923693044;% 12 NPY
    2.83370985669454e-05;% 13 NPY vsm
    0.290000000000000; %vol a
    0.440000000000000; %vol c
    0.270000000000000; %vol v
    1;1;1; %flows
    2.11708250620889;% n a
    2.28833185479117;%n c
    1.47095919879372;%n v
    1.13809920000000;%n t
    0.3; %Cn
    0.3; %Cg
    0.3; %Cc
    0];%ATP
%options.x0 = sol.x(end,:).';
TE = 20*10^-3;       B0 = 7;
Constants1 = [7.35304240575339e-06, 0.000581299567725387, 2.83400100550909e-05,...
    Ca_start, 0, Con, HbO_0, HbR_0, SaO2_0, ScO2_0, SvO2_0, TE, B0,...
    cCO2in, 0.35];
Constants2 = [7.35304240575339e-06, 0.000581299567725387, 2.83400100550909e-05,...
    Ca_start, 20, Con, HbO_0, HbR_0, SaO2_0, ScO2_0, SvO2_0, TE, B0,...
    cCO2in, 0.35];
%%
% alter simulation tolerances, DAE solver can not handle the default values
options.atol = 1e-6;
options.rtol = 1e-12;
solReal1 = simulate_Dyn_H22_A3r(0:1:200,theta,Constants1,[],options);
options.x0 = solReal1.x(end,:).';
solReal2 = simulate_Dyn_H22_A3r(0:.1:50,theta,Constants2,[],options);

figure(5)
plot(solReal2.t,solReal2.y(:,2),'r')







