clc;clear;close all;
% Plot Setup Variables
ExpNum=1;
PlotSimFlag=1;
plotData='';
% Plotting Code - No Need To Change
PWD=pwd;
addpath(strcat(pwd,'\Other_Scripts'));
cd(strcat(PWD,'\Exp_Data'))
Data=load(strcat('Exp',num2str(ExpNum),'.mat'));
Data=Data.Data;
cd(PWD);
titles=Data{2,end};
[rows,cols]=size(Data);
maxcbf=0;
mincbf=5;
maxbold=0;
minbold=5;
for i=1:cols-1
    maxcbf=max([maxcbf max(Data{1,i}+Data{2,i})]);
    mincbf=min([mincbf min(Data{1,i}-Data{2,i})]);
    maxbold=max([maxbold max(Data{4,i}+Data{5,i})]);
    minbold=min([minbold min(Data{4,i}-Data{5,i})]);
end
bounds=[(mincbf-0.05*maxcbf) (maxcbf+0.05*maxcbf);(minbold-0.05*maxbold) (maxbold+0.05*maxbold)];
figure(ExpNum)
for i=1:2
    for j=1:(cols-1)
        subplot(2,(cols-1),(cols-1)*(i-1)+j)
        errorbar(Data{1,end},Data{1+3*(i-1),j},Data{1+3*(i-1)+1,j},'k*')
        hold on
        title(titles{j},'FontSize',12,'FontWeight','Bold')
        xlabel('time (s)','FontSize',12,'FontWeight','Bold')
        if i==1
            ylabel('cbf','FontSize',12,'FontWeight','Bold')
        else
            ylabel('%\DeltaBOLD','FontSize',12,'FontWeight','Bold')
        end
        ylim([bounds(i,1) bounds(i,2)]);
        xlim([0 140])
    end
end

%%
if PlotSimFlag==1
    [~, ~, Con, stimend, ~] = optsetupfunction(1);
    cd(strcat(PWD,'\Optimized_Params'))
    if isempty(plotData)==1
        [file,path] = uigetfile;
        cd(path)
        load(file,'X');
    else
        load(plotData,'X');
    end
    cd(PWD)
    theta=X';
    options = amioption('sensi',0,...
        'maxsteps',1e4);
    options.sensi = 0;
    Ca_start = 10;
    sol = simulate_SSmodel(inf,theta(4:end),[Ca_start,Con],[],options);
    ssArt = sol.y(1);
    HbO_0 = sol.y(2);
    HbR_0 = sol.y(3);
    SaO2_0 = sol.y(4);
    ScO2_0 = sol.y(5);
    SvO2_0 = sol.y(6);
    options.x0 = sol.x(end,:).';
    TE = 20*10^-3;       B0 = 7;
    Constants = [sol.x(end,[11 9 13]), Ca_start, stimend(3)-stimend(2), Con, HbO_0, HbR_0, SaO2_0, ScO2_0, SvO2_0, TE, B0];
    options.atol = 1e-6;
    options.rtol = 1e-12;
    solReal1 = simulate_H22(Data{1,cols}(1:12),theta, Constants, [], options);
    solReal2 = simulate_H22(0:0.1:100,theta, Constants, [], options);
    
    figure(1)
    subplot(2,(cols-1),4)
    plot(solReal1.t+30,solReal1.y(:,6),'r')
    R1=plot(solReal2.t+30,solReal2.y(:,6),'r');
    R1.Color = [R1.Color 0.4];
    subplot(2,(cols-1),1)
    plot(solReal1.t+30,solReal1.y(:,5),'r')    
    R1=plot(solReal2.t+30,solReal2.y(:,5),'r');
    R1.Color = [R1.Color 0.4];
end