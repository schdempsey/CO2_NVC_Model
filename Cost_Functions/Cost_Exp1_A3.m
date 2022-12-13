function [f,c,gf,gc] = Cost_Exp1_A3(theta,Data,Con,stimend,FID)
    logL=0;
    dlogL = zeros(numel(theta),1);
    options = amioption('sensi',0,...
        'maxsteps',1e4);
    options.sensi = 0;
    options.nmaxevent = 1e3;
    options.atol = 1e-6;
    options.rtol = 1e-12;
    stimTime=stimend(3)-stimend(2); %stimtime
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
    C5=cCO2in*1.05;
    %C5=0.5;
    Constants2(5)=0; %No neuronal stim, but CO2 can change G because Tog=1 
    Constants2(20)=C5; %new CO2 concentration
    solSSC5 = simulate_Dyn_H22_A3(inf,theta,Constants2,[],options); %sim no stim with new C inf 
    %solSSC5r = simulate_Dyn_H22_A3(0:.1:110,theta,Constants2,[],options); %sim no stim with new C growth
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
    %%
    if (solC0.status<0) && (solC5.status<0)
        logL = NaN;
    else
        costCBFC0 = sum((solC0.y(1:100:1101,2) - Data{1,1}(4:end)).^2./(Data{3,1}(4:end).^2));
        costBOLDC0 = sum((solC0.y(1:100:1101,3) - Data{4,1}(4:end)).^2./(Data{6,1}(4:end).^2));
        costCBFC5 = sum((solC5.y(1:100:1101,2) - Data{1,2}(4:end)).^2./(Data{3,2}(4:end).^2));
        costBOLDC5 = sum((solC5.y(1:100:1101,3) - Data{4,2}(4:end)).^2./(Data{6,2}(4:end).^2));
        logL = logL + costCBFC0 + costBOLDC0 + costCBFC5 + costBOLDC5;
    end
    
    f = logL;
    gf = dlogL;
    c = [];
    gc = [];
    %% MCMC related, save parameters to file
    if nargin == 5 && logL < chi2inv(0.95,11) 
        fprintf(FID,'%4.10f %10.10f ',[f, theta']); fprintf(FID,'\n');
    end
end