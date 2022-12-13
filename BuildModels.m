addpath('Model_Files')
addpath('Model_syms')
path = [pwd, '/Model_Files/'];

%amiwrap('H22','H22_syms',path);
%amiwrap('SSM_H22_A3','SSM_H22_A3_syms',path)
amiwrap('Dyn_H22_A3','Dyn_H22_A3_syms',path)