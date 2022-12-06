%% Constatants needed for the circuit part of the model
function [Constants] = GenerateConstants()

pO2_01 = 81.2;      pO2_12 = 59.7;      pO2_23 = 39.6;      pO2_34 = 41.3;
pO2_t = 22.4;       pO2_femart = 85.6;

[g1,g2,g3,gs,CMRO2_0,CO2_l] = BarretPO2constants(pO2_01, pO2_12, pO2_23, pO2_34, pO2_t, pO2_femart);
Constants = [g1,g2,g3,gs,CMRO2_0,CO2_l,pO2_femart];

end