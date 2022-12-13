clear;clc;
addpath(strcat(pwd,'\Model_Files'));
addpath(strcat(pwd,'\Other_Scripts'));
options = amioption('sensi',0,...
    'maxsteps',1e4);
options.sensi = 0;
options.nmaxevent = 1e3;
options.atol = 1e-6;
options.rtol = 1e-12;
%% Input Parameters
Constants=load('Constants.mat');
Con=Constants.Constants;
theta=[2.11020855546808;1.15116690047892;-0.998721418427387;3;1.45501023492916;-4.46362660461291;-0.724077891135583;-3.99573240520881;0.523924736331482;2.14718556419683;4.29758639199676;0.111300152345884;1.67436525595713;3.98614448351813;0.121508516284357;0.363416836400271;4.37969324599745;4.49641225986757;-0.357181797944409;2.19545317336248;-4.05935598126191;4.50000000000000;-1.84127655079188;-3.74850175777620;-11.6970919520096;-4.50000000000000;-4.50000000000000;-0.721186153425439;4.48894055158979;0.292103231526475;0.122718614591054;0.161232435616766;1;2.38458129653135;2.59211603267482;-0.698430359222535;-4.09909370339094;-0.830776016508415;7.36591832140683;2.06101156570888;3.16106323796090;3.85024600387078;0.944743801332074;1.71805062728506;0.0276775850927731];
stimTime=60;
%% Steady state simulation
Ca_start = 10; %Constant Calcium Production
cCO2in=0.3; %percent CO2 in
TE = 20*10^-3; %Time of Echo
B0 = 7; %Magnetic Field Strength
Tog=0; %No G
Constants1 = [0,0,0,Ca_start,0,Con,0,0,0,0,0,TE,B0,cCO2in,0.35,Tog,0.29,0.44,0.27];
solSSCO = simulate_Dyn_H22_A3(inf,theta,Constants1,[],options); %Run Initial Steady State
%%  Get new C5 SS inputs
HbO_0 = solSSCO.y(6);
HbR_0 = solSSCO.y(7);
SaO2_0 = solSSCO.y(8);
ScO2_0 = solSSCO.y(9);
SvO2_0 = solSSCO.y(10);
cCO2n0 = solSSCO.y(11);
Tog=1; %Now activate G
Constants2 = [solSSCO.x(11),solSSCO.x(9),solSSCO.x(13),Ca_start,stimTime,Con,HbO_0,HbR_0,SaO2_0,ScO2_0,SvO2_0,TE,B0,cCO2in,cCO2n0,Tog,solSSCO.x(14),solSSCO.x(15),solSSCO.x(16)];
options.x0 = solSSCO.x(end,:).'; %add steady state as initial conditions
%% Normocapnia Stimulation
solC0 = simulate_Dyn_H22_A3((0:.1:110)',theta,Constants2,[],options); %sim stim with C0 
%% New CO2 Steady State Sim- 5% increase
%C5=cCO2in*1.05;
C5=0.5;
Constants2(5)=0; %No neuronal stim, but CO2 can change G because Tog=1 
Constants2(20)=C5; %new CO2 concentration
solSSC5 = simulate_Dyn_H22_A3(inf,theta,Constants2,[],options); %sim no stim with new C inf 
solSSC5r = simulate_Dyn_H22_A3(0:.1:110,theta,Constants2,[],options); %sim no stim with new C growth
%% Get new C5 SS BOLD ONLY inputs
HbO_0 = solSSC5.y(6);
HbR_0 = solSSC5.y(7);
SaO2_0 = solSSC5.y(8);
ScO2_0 = solSSC5.y(9);
SvO2_0 = solSSC5.y(10);
V1=solSSC5.x(14);
V2=solSSC5.x(15);
V3=solSSC5.x(16);
Constants3=Constants2;
Constants3(5)=stimTime;
Constants3(13:17) = [HbO_0,HbR_0,SaO2_0,ScO2_0,SvO2_0];
Constants3(23:25) = [V1,V2,V3];
%% New CO2 Steady State Sim
options.x0 = solSSC5.x(end,:).'; %input 5% CO2 SS
solC5 = simulate_Dyn_H22_A3((0:.1:110)',theta,Constants3,[],options);


figure(1)
subplot(2,3,1)
plot(solC0.t(1:100:1101),solC0.y(1:100:1101,2),'r')
subplot(2,3,4)
plot(solC0.t(1:100:1101),solC0.y(1:100:1101,3),'r')
hold on
subplot(2,3,2)
plot(solSSC5r.t(1:100:1101),solSSC5r.y(1:100:1101,2),'r')
subplot(2,3,5)
plot(solSSC5r.t(1:100:1101),solSSC5r.y(1:100:1101,3),'r')
hold on
subplot(2,3,3)
plot(solC5.t(1:100:1101),solC5.y(1:100:1101,2),'r') 
subplot(2,3,6)
plot(solC5.t(1:100:1101),solC5.y(1:100:1101,3),'r')







