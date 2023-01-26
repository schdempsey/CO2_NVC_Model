% help function [objectiveFunction, Data, Constants, stimparams, X] = optsetupfunction(choice)
% Outputs:
%   objectiveFunction - callable pointer to the cost-function, with all
%   additional input arguments specified
function [objectiveFunction, Data, Constants, stimparams, X] = optsetupfunction(choice)
    Constants=load('Constants.mat');
    Constants=Constants.Constants;
    stimparams=60; %stimulation time
    PWD=pwd;
    cd(strcat(PWD,'\Exp_Data'))
    Data=load(strcat('Exp',num2str(choice),'.mat'));
    Data=Data.Data;
    cd(PWD);
    %Eventually add cases for each model as well, A3,A2,A1 but for now all
    %A3.
    switch choice
        case 1
            X = zeros(45,1);
            objectiveFunction=@(X) Cost_Exp1_A3(X,Data,Constants,stimparams);
        case 2
            X = zeros(45,1);
            objectiveFunction=@(X) Cost_Exp2_A3(X,Data,Constants,stimparams);
        case 3
            X = zeros(45,1);
            objectiveFunction=@(X) Cost_Exp3_A3(X,Data,Constants,stimparams);
    end
end