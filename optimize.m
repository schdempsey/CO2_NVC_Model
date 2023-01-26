close all;clear;
PWD=pwd;
addpath(strcat(pwd,'\Cost_Functions'));
addpath(strcat(pwd,'\Exp_Data'));
addpath(strcat(pwd,'\Model_Files'));
addpath(strcat(pwd,'\Other_Scripts'));
% Script to perform an optimization using ESS-opt
%% Load inital data structures
ExpNum=1;
Model='A3';
[objectiveFunction, Data, Constants, stimend, X] = optsetupfunction(ExpNum);
%% Generate parameter bounds
lb = -4.5*ones(size(X));
ub =4.5*ones(size(X));
ub([4 5 6 7 8 9]) = 3;                       %KPF & KPINF
lb([25 37]) = -12;                           %
lb([10 11 12 13 14 15]) = log10(1/0.78);     %ksink 	    
ub(1:3) = 2.5;  
lb(22) = 0;                                  %sinkNO
lb(30:32) = [0.01 log10(1.3) log10(1.4)];    %K123 circuit
lb(33:35) = [1 1 1];                         %vis123 circuit
ub(30:32) = [log10(2) log10(2) 8];           %K123 circuit
ub(33:35) = [2 2.7 2.7];  %vis123 circuit

ub(37:45) = 8; %new Param upper bounds
lb(37:45) = -8; %new Param lower bounds

problem.x_L       = lb; % essOPT uses a problem structure where crucial information is specified
problem.x_U       = ub;
problem.vtr=-100;
if sum(X)==0
    X=(lb+ub)/2;
end
problem.x_0=X; %Initial params set from either optimized data or scratch.
%% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse   =100;       %100; %500; %5; %
opts.maxtime    = 2000;%200;      % MAX-Time of optmization, i.e how long the optimization will last
opts.maxeval    = 1e8;      % max number of evals, i.e cost function calls
opts.log_var    = [];    
opts.local.solver = 'dhc';%dhc'; %'fmincon'; %'nl2sol'; %'mix'; 
opts.local.finish = opts.local.solver; %uses the local solver to check the best p-vector
opts.local.bestx = 0;      
opts.local.balance = 2;   %how far from startguess the local search will push the params, 0.5 default
opts.local.n1   = 5;%2;        %Number of iterations before applying local search for the 1st time (Default 1)
opts.local.n2   = 5;%1;        %Minimum number of iterations in the global phase between 2 local calls (Default 10) 
problem.f       = 'meigoDummy'; % calls function that sets up the cost function call
%% MEIGO OPTIONS II (FOR ESS AND MULTISTART):
opts.local.iterprint = 1; % prints what going on during optimization
%% MEIGO OPTIONS III (FOR ESS ONLY):
opts.dim_refset   = 10;
%% OPTIONS AUTOMATICALLY SET AS A RESULT OF PREVIOUS OPTIONS:
if(strcmp(opts.local.solver,'fmincon'))
    opts.local.use_gradient_for_finish = 1; %DW: provide gradient to fmincon
else
    opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
end
opts.local.check_gradient_for_finish = 0; %DW: gradient checker
%% Solve
warning('off','all') 
optim_algorithm = 'ess'; % 'multistart'; %  'cess'; 
Results = MEIGO(problem,opts,optim_algorithm,objectiveFunction); % Run the optimization
%% Save results
fitting_speed     = Results.time(end);
best_fs           = Results.fbest;
parameters_ess    = Results.xbest';
% best parameter vector stored in X
X=parameters_ess;
w = warning ('on','all'); 

%Save the results with date tag, model tag, and run number for that day
tag1=datetime("today");
[y,m,d] = ymd(tag1);
tag1=strcat(num2str(y),num2str(m,'%02.f'),num2str(d,'%02.f'));
tag2=Model;
filename=strcat(tag1,'_','Exp',num2str(ExpNum),tag2,'_run');
%filename=strcat(tag1,'_','Exp',num2str(1),tag2,'_run');
cd(strcat(PWD,'\Optimized_Params\Exp',num2str(ExpNum)))
%cd(strcat(PWD,'\Optimized_Params\Exp',num2str(1)))
DIR=dir;
runNum=0;
for j=3:(length(DIR))
    if strcmp(DIR(j).name(1:19),filename)==1
        runNum=max([runNum str2num(DIR(j).name(20:21))]);
    end
end
runNum=runNum+1;
filename=strcat(filename,num2str(runNum,'%02.f'),'.mat');
save(filename,'X')
cd(PWD)
