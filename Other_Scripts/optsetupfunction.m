% help function [objectiveFunction, Data, Constants, stimparams, X] = optsetupfunction(choice)
% Outputs:
%   objectiveFunction - callable pointer to the cost-function, with all
%   additional input arguments specified
function [objectiveFunction, Data, Constants, stimparams, X] = optsetupfunction(choice)
    Constants=load('Constants.mat');
    Constants=Constants.Constants;
    stimparams=[0 30 90 14 10];
    PWD=pwd;
    cd(strcat(PWD,'\Exp_Data'))
    Data=load(strcat('Exp',num2str(choice),'.mat'));
    Data=Data.Data;
    cd(PWD);
    switch choice
        case 1
            X = zeros(37,1);
            objectiveFunction=@(X) Cost_Exp1(X,Data,Constants,stimparams);
        case 2
            X = zeros(37,1);
            objectiveFunction=@(X) Cost_Exp2(X,Data,Constants,stimparams);
        case 3
            X = ones(37,1);
            objectiveFunction=@(X) Cost_Exp3(X,Data,Constants,stimparams);
    end
end