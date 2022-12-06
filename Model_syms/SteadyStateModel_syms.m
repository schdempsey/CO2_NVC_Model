function [model] = SteadyStateModel_syms()

% set the parametrisation of the problem options are 'log', 'log10' and

model.param = 'log10';

%%
% STATES
% create state syms 
syms N_NO N_NPY N_Pyr Ca_NO Ca_NPY Ca_Pyr AA PGE2 PGE2vsm NO NOvsm NPY NPYvsm V1 V2 V3 f1 f2 f3 nO2_1 nO2_2 nO2_3 nO2_t 
% create state vector
model.sym.x = [N_NO, N_NPY, N_Pyr, Ca_NO, Ca_NPY, Ca_Pyr, AA, PGE2, PGE2vsm, NO, NOvsm, NPY, NPYvsm, V1, V2, V3, f1, f2, f3, nO2_1, nO2_2, nO2_3, nO2_t];
%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms                                                                
syms kPF1 kPF2 kIN kIN2 kINF kINF2 sinkN_NO sinkN_NPY sinkN_Pyr sinkCa_NO sinkCa_NPY sinkCa_Pyr kPL kCOX kPGE2 sinkPGE2 kNOS kNO sinkNO kNPY Vmax Km sinkNPY ky1 ky2 ky3  K1 K2 K3 vis1 vis2 vis3 kscalemet Km2
% create parameter vector 
model.sym.p = [kPF1, kPF2, kIN, kIN2, kINF, kINF2, sinkN_NO, sinkN_NPY, sinkN_Pyr, sinkCa_NO, sinkCa_NPY, sinkCa_Pyr, kPL, kCOX, kPGE2, sinkPGE2, kNOS, kNO, sinkNO, kNPY, Vmax, Km, sinkNPY, ky1, ky2, ky3, K1, K2, K3, vis1, vis2, vis3, kscalemet, Km2];
%%  
% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
syms kCa g_1 g_2 g_3 g_s CMRO2_0 CO2_l pO2_femart
% create parameter vector 
model.sym.k = [kCa, g_1, g_2, g_3, g_s, CMRO2_0, CO2_l, pO2_femart];
%%
% SYSTEM EQUATIONS
% create symbolic variable for time
syms t 
model.sym.xdot = sym(zeros(size(model.sym.x)));

%% Neuronal dynamics (zero in steady-state)
model.sym.xdot(1) = 0;
model.sym.xdot(2) = 0;
model.sym.xdot(3) = 0;

%% Calcium dynamics
model.sym.xdot(4) = kCa-sinkCa_NO*Ca_NO;
model.sym.xdot(5) = kCa-sinkCa_NPY*Ca_NPY;
model.sym.xdot(6) = kCa-sinkCa_Pyr*Ca_Pyr;

%% Pyramidal intracellular signaling AA->PGE2
model.sym.xdot(7) = kPL*Ca_Pyr-kCOX*AA/(Km2+AA);
model.sym.xdot(8) = kCOX*AA/(Km2+AA)-kPGE2*PGE2;
model.sym.xdot(9) = kPGE2*PGE2-sinkPGE2*PGE2vsm;

%% GABAergic NO signaling
model.sym.xdot(10) = kNOS*Ca_NO-kNO*NO;
model.sym.xdot(11) = kNO*NO-sinkNO*NOvsm;

%% GABAergic NPY signaling
model.sym.xdot(12) = kNPY*Ca_NPY-Vmax*NPY/(Km+NPY);
model.sym.xdot(13) = Vmax*NPY/(Km+NPY)-sinkNPY*NPYvsm;

Artdia = ky1*NOvsm +ky2*PGE2vsm -ky3*NPYvsm;
%% Circuit
model.sym.xdot(14) = 0;
model.sym.xdot(15) = 0;
model.sym.xdot(16) = 0;

f0 = 1;
model.sym.xdot(17) = 0;
model.sym.xdot(18) = 0;
model.sym.xdot(19) = 0;

% CBF = (V1*(f0 + f1) + V2*(f1 + f2) + V3*(f2 + f3))/(2*(V1 + V2 + V3));

%% Barret Hb and pressure (mM deafault unit??)
CO2_max = 9.26;           
V_t = 34.8;
p50 = 36;           h = 2.6;
sigma_O2 = 1.46*10^-3;

CO2_0 = CO2_max/(10^(log10(pO2_femart/p50)/(-1/h))+1);
CO2_01 = CO2_0 - CO2_l;
CO2_12 = nO2_1/V1;
CO2_23 = nO2_2/V2;
CO2_34 = nO2_3/V3;
CO2_t = nO2_t/V_t;        

pO2_01 = p50*((CO2_max/CO2_01) -1)^(-1/h);   
pO2_12 = p50*((CO2_max/CO2_12) -1)^(-1/h);
pO2_23 = p50*((CO2_max/CO2_23) -1)^(-1/h);
pO2_34 = p50*((CO2_max/CO2_34) -1)^(-1/h);

pO2_1 = (pO2_01 + pO2_12)/2;
pO2_2 = (pO2_12 + pO2_23)/2;
pO2_3 = (pO2_23 + pO2_34)/2;
pO2_t = CO2_t/sigma_O2;

jO2_1 = g_1*(pO2_1 - pO2_t);
jO2_2 = g_2*(pO2_2 - pO2_t);
jO2_3 = g_3*(pO2_3 - pO2_t);
jO2_s = g_s*(pO2_1 - pO2_3);

model.sym.xdot(20) = f0*CO2_01 - f1*CO2_12 - jO2_1 - jO2_s;
model.sym.xdot(21) = f1*CO2_12 - f2*CO2_23 - jO2_2;
model.sym.xdot(22) = f2*CO2_23 - f3*CO2_34 - jO2_3 + jO2_s;
model.sym.xdot(23) = jO2_1 + jO2_2 + jO2_3 - CMRO2_0;

SaO2 = ((CO2_01 + CO2_12)/2)/(CO2_max);
ScO2 = ((CO2_12 + CO2_23)/2)/(CO2_max);
SvO2 = ((CO2_23 + CO2_34)/2)/(CO2_max);

HbOa = V1*SaO2;     HbRa = V1*(1-SaO2);          
HbOc = V2*ScO2;     HbRc = V2*(1-ScO2);           
HbOv = V3*SvO2;     HbRv = V3*(1-SvO2);        

HbO = HbOa + HbOc + HbOv;
HbR = HbRa + HbRc + HbRv;  
%% INITIAL CONDITIONS
model.sym.x0(1:13) = 0;
model.sym.x0(14:19) = [0.29; 0.44; 0.27; 1; 1; 1];
model.sym.x0(20:23) = 1; %[2.1171; 2.7284; 1.471; 1.1381];

model.sym.dx0 = sym(zeros(size(model.sym.x)));

% OBSERVALES
model.sym.y = sym(zeros(6,1));
model.sym.y(1) = Artdia;
model.sym.y(2) = HbO;
model.sym.y(3) = HbR;
model.sym.y(4) = SaO2;
model.sym.y(5) = ScO2;
model.sym.y(6) = SvO2;
