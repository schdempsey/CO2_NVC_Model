function [f,c,gf,gc] = Cost_Exp1(theta,Data,Con,stimend,FID)
    logL=0;
    dlogL = zeros(numel(theta),1);
    options = amioption('sensi',0,...
        'maxsteps',1e4);
    options.sensi = 0;
    options.nmaxevent = 1e3;
    %% Steady state simulation
    Ca_start = 10;
    % steady state simulation
    sol = simulate_SSmodel(inf,theta(4:end),[Ca_start,Con],[],options);
    % assaign values to constants in the stimulation simulation
    HbO_0 = sol.y(2);
    HbR_0 = sol.y(3);
    SaO2_0 = sol.y(4);
    ScO2_0 = sol.y(5);
    SvO2_0 = sol.y(6);
    options.x0 = sol.x(end,:).';
    TE = 20*10^-3;       B0 = 7;
    Constants = [sol.x(end,[11 9 13]), Ca_start, stimend(3)-stimend(2), Con, HbO_0, HbR_0, SaO2_0, ScO2_0, SvO2_0, TE, B0];

    % alter simulation tolerances, DAE solver can not handle the default values
    options.atol = 1e-6;
    options.rtol = 1e-12;

    %% Simulations
        solReal = simulate_H22((0:10:110)',theta, Constants, [], options);
        options.atol = 1e-5;
        options.rtol = 1e-8;
        if (solReal.status<0)
            logL = NaN;
        else
            costCBF = sum((solReal.y(:,5) - Data{1,1}(4:end)).^2./(Data{3,1}(4:end).^2));
            costBOLD = sum((solReal.y(:,6) - Data{4,1}(4:end)).^2./(Data{6,1}(4:end).^2));
            logL = logL + costCBF+costBOLD;
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