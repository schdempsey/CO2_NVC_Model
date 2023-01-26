clc;clear;%close all;
%%%%%%%%%%%%%%%%%%%%%%
% Plot Setup Variables
ExpNum=1;
PlotSimFlag=1; %Make zero if you want to just plot the experimental data
plotMetabolites=1;
datafilename=''; %If you want to load a specific filename you can put it 
    %here, otherwise you can pick from the file browser if left blank
%%%%%%%%%%%%%%%%%%%%%%
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
boundLim=0.05; %same for now, but could be diff for CBF and BOLD
bounds=[(mincbf-boundLim*maxcbf) (maxcbf+boundLim*maxcbf);...
    (minbold-boundLim*maxbold) (maxbold+boundLim*maxbold)];
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
        %xlim([0 140])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Sim Results and Plot on Top of Experimental REsults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PlotSimFlag==1
    [~, ~, Con, stimend, ~] = optsetupfunction(ExpNum);
    cd(strcat(PWD,'\Optimized_Params'))
    if isempty(datafilename)==1
        [file,path] = uigetfile;
        cd(path)
        load(file,'X');
    else
        load(datafilename,'X');
    end
    cd(PWD)
    theta=X'; %optimized params
    options = amioption('sensi',0,...
        'maxsteps',1e6);
    options.sensi = 0;
    options.nmaxevent = 1e3;
    options.atol = 1e-6;
    options.rtol = 1e-12;
    stimTime=60; %stimtime, do not change.
    simLen=300;%110; %default 110, can make longer if you want to see more
    fontS=10; %font size for all figs.
    %% Steady state simulation
    Ca_start = 10; %Constant Calcium Production
    cCO2in=1; % CO2 in normalized to 1
    C5=cCO2in*1.05; %increase of 5%. NOTE: This is only possible because CO2 partial pressures are linear with CO2 amount
    C10=cCO2in*1.1; %increase of 10%
    TE = 10*10^-3; %Time of Echo 10 for H22
    B0 = 9.4; %Magnetic Field Strength 9.4T for H22
    Tog=0; %No G activation at steady state. Will just initialize the oxygen and neuron params
    Constants1 = [0,0,0,Ca_start,0,Con,0,0,0,0,0,TE,B0,cCO2in,0.35,Tog,0.29,0.44,0.27];
    solSSCO = simulate_Dyn_H22_A3(inf,theta,Constants1,[],options); %Run Initial Steady State
    %%  Get new C5 SS inputs
    HbO_0 = solSSCO.y(3);
    HbR_0 = solSSCO.y(4);
    SaO2_0 = solSSCO.y(5);
    ScO2_0 = solSSCO.y(6);
    SvO2_0 = solSSCO.y(7);
    cCO2n0 = solSSCO.y(8);
    Tog=1; %Now activate G
    Constants2 = [solSSCO.x(11),solSSCO.x(9),solSSCO.x(13),Ca_start,...
        stimTime,Con,HbO_0,HbR_0,SaO2_0,ScO2_0,SvO2_0,TE,B0,cCO2in,...
        cCO2n0,Tog,solSSCO.x(14),solSSCO.x(15),solSSCO.x(16)];
    options.x0 = solSSCO.x(end,:).'; %Make steady state initial conditions
    %% Normocapnia Stimulation
    solC0 = simulate_Dyn_H22_A3((0:.1:simLen)',theta,Constants2,[],options); %sim stim with C0 
    %% New CO2 Steady State Sim at 5% increase
    Constants2(5)=0; %No neuronal stim, but CO2 can change G and dilation because Tog=1 
    Constants2(20)=C5; %new CO2 concentration
    solSSC5 = simulate_Dyn_H22_A3(0:.1:500,theta,Constants2,[],options); %sim no stim with new C inf 
    %% Get new C5 SS BOLD ONLY inputs
    HbO_0 = solSSC5.y(end,3);
    HbR_0 = solSSC5.y(end,4);
    SaO2_0 = solSSC5.y(end,5);
    ScO2_0 = solSSC5.y(end,6);
    SvO2_0 = solSSC5.y(end,7);
    V1=solSSC5.x(end,14);
    V2=solSSC5.x(end,15);
    V3=solSSC5.x(end,16);
    Constants3=Constants2;
    Constants3(5)=stimTime;
    Constants3(13:17) = [HbO_0,HbR_0,SaO2_0,ScO2_0,SvO2_0];
    Constants3(23:25) = [V1,V2,V3];
    %% New CO2 Steady State Sim
    options.x0 = solSSC5.x(end,:).'; %input 5% CO2 SS
    solC5 = simulate_Dyn_H22_A3((0:.1:simLen)',theta,Constants3,[],options);
    %% New CO2 Steady State Sim- 10% increase PREDICTION
    Constants3(5)=0; %No neuronal stim, but CO2 can change G because Tog=1 
    Constants3(20)=C10; %new CO2 concentration
    %solSSC10 = simulate_Dyn_H22_A3(inf,theta,Constants3,[],options); %sim no stim with new C inf 
    solSSC10 = simulate_Dyn_H22_A3(0:.1:5000,theta,Constants3,[],options); 
    %sim no stim with new C for 5000 seconds because the inf option might have failed
    %% Get new C5 SS BOLD ONLY inputs
    HbO_0 = solSSC10.y(end,3);
    HbR_0 = solSSC10.y(end,4);
    SaO2_0 = solSSC10.y(end,5);
    ScO2_0 = solSSC10.y(end,6);
    SvO2_0 = solSSC10.y(end,7);
    V1=solSSC10.x(end,14);
    V2=solSSC10.x(end,15);
    V3=solSSC10.x(end,16);
    Constants4=Constants3;
    Constants4(5)=stimTime;
    Constants4(13:17) = [HbO_0,HbR_0,SaO2_0,ScO2_0,SvO2_0];
    Constants4(23:25) = [V1,V2,V3];
    %% New CO2 Steady State Sim
    options.x0 = solSSC10.x(end,:).'; %input 5% CO2 SS
    solC10 = simulate_Dyn_H22_A3((0:.1:simLen)',theta,Constants4,[],options);
    %% Plot Everything
    step=100;
    figure(1)
    subplot(2,3,1)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),1),'r')
    subplot(2,3,4)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),2),'r')
    subplot(2,3,2)
    plot(solC5.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),1),'r')
    subplot(2,3,5)
    plot(solC5.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),2),'r')
    subplot(2,3,3)
    plot(solC10.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),1),'b')
    subplot(2,3,6)
    plot(solC10.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),2),'b')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%METABOLITES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotMetabolites==1
    %% Make a Figure to Analyze Intermediates for Discussion and learning
    figure(2)
    step=5; %time step for data
    rows=6; %rows to add to massive subplot.
    %% Vasoactive Effect
    subplot(rows,3,1)
    plot(solC0.t(1:step:(simLen*10+1))+30,10^(theta(27))*(solC0.x(1:step:(simLen*10+1),11)-solC0.x(1,11)),'Color',[0.8500 0.3250 0.0980]) %NOvsm
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,10^(theta(28))*(solC0.x(1:step:(simLen*10+1),9)-solC0.x(1,9)),'Color',[0.4940 0.1840 0.5560]); %PGE vsm
    plot(solC0.t(1:step:(simLen*10+1))+30,-10^(theta(29))*(solC0.x(1:step:(simLen*10+1),13)-solC0.x(1,13)),'Color',[0.3010 0.7450 0.9330]) %NPY vsm
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),12),'k--','LineWidth',1.5) %ArtDia
    a=legend('NOvsm','PGEvsm','NPYvsm','Net Activity Function');
    set(a,'FontSize',fontS)
    ylabel('Dilation fn Normocapnia','FontSize',fontS)
    %ylim([-0.1 1.6])
    subplot(rows,3,2)
    plot(solC0.t(1:step:(simLen*10+1))+30,10^(theta(27))*(solC5.x(1:step:(simLen*10+1),11)-solC0.x(1,11)),'Color',[0.8500 0.3250 0.0980])
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,10^(theta(28))*(solC5.x(1:step:(simLen*10+1),9)-solC0.x(1,9)),'Color',[0.4940 0.1840 0.5560]);
    plot(solC0.t(1:step:(simLen*10+1))+30,-10^(theta(29))*(solC5.x(1:step:(simLen*10+1),13)-solC0.x(1,13)),'Color',[0.3010 0.7450 0.9330])
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),12),'k--','LineWidth',1.5)
    ylabel('Dilation fn 5%CO_2','FontSize',fontS)
    %ylim([-0.1 1.6])
    subplot(rows,3,3)
    plot(solC0.t(1:step:(simLen*10+1))+30,10^(theta(27))*(solC10.x(1:step:(simLen*10+1),11)-solC0.x(1,11)),'Color',[0.8500 0.3250 0.0980])
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,10^(theta(28))*(solC10.x(1:step:(simLen*10+1),9)-solC0.x(1,9)),'Color',[0.4940 0.1840 0.5560]);
    plot(solC0.t(1:step:(simLen*10+1))+30,-10^(theta(29))*(solC10.x(1:step:(simLen*10+1),13)-solC0.x(1,13)),'Color',[0.3010 0.7450 0.9330])
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),12),'k--','LineWidth',1.5) 
    ylabel('Dilation fn 10%CO_2','FontSize',fontS)
    %ylim([-0.1 1.6])
     %% Volume Response
    subplot(rows,3,4)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),14),'r') %V art
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),14),'g') %V ven
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),14),'b') %V cap
    a=legend('Normo','5% CO_2','10% CO_2');
    set(a,'FontSize',fontS)
    ylabel('v_a','FontSize',fontS);
    ylim([0.2 0.6])
    subplot(rows,3,5)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),15),'r') 
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),15),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),15),'b') 
    ylabel('v_c','FontSize',fontS);
    ylim([0.2 0.6])
    subplot(rows,3,6)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),16),'r') 
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),16),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),16),'b')
    ylabel('v_v','FontSize',fontS);
    ylim([0.2 0.6])
    %% CO2 Concentration Response
    subplot(rows,3,7)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),24),'r') %CO_2 Neuron
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),24),'g') %CO_2 AStro
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),24),'b') %CO_2 Cap
    ylabel('[CO_2]_n','FontSize',fontS);
    %ylim([0.8 2.7])
    subplot(rows,3,8)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),25),'r') 
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),25),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),25),'b')
    ylabel('[CO_2]_a','FontSize',fontS);
    %ylim([0.8 2.7])
    subplot(rows,3,9)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),26)./solC0.x(1:step:(simLen*10+1),15),'r')
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),26)./solC5.x(1:step:(simLen*10+1),15),'g')
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),26)./solC10.x(1:step:(simLen*10+1),15),'b')
    ylabel('[CO_2]_c','FontSize',fontS);
    %ylim([0.8 2.7])
    %% CMRO2 and Knockout Response
    subplot(rows,3,10)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),11),'r') %CMRO2
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),11),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),11),'b') 
    ylabel('CMRO_2','FontSize',fontS);
    %ylim([0.8 2.7])
    subplot(rows,3,11)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),9),'r') %KO NO
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),9),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),9),'b') 
    ylabel('NO KO','FontSize',fontS);
    %ylim([0.8 2.7])
    subplot(rows,3,12)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),10),'r') %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),10),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),10),'b') 
    ylabel('NPY KO','FontSize',fontS)
    %ylim([0.8 2.7])
    %% ATP 
    subplot(rows,3,13)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),27),'r') %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),27),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),27),'b') 
    ylabel('[ATP]','FontSize',fontS)
    subplot(rows,3,14)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),16),'r') %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),16),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),16),'b') 
    ylabel('pO_2 tissue','FontSize',fontS)
    subplot(rows,3,16)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),13),'Color',[0.4660 0.6740 0.1880]) %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),14),'Color',[0.6350 0.0780 0.1840]) 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.y(1:step:(simLen*10+1),15),'Color',[0 0.4470 0.7410]) 
    a=legend('HbO_T','HbO_O','HbO_R');
    set(a,'FontSize',fontS);
    ylabel('HbO Normo','FontSize',fontS)
    subplot(rows,3,17)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),13),'Color',[0.4660 0.6740 0.1880]) %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),14),'Color',[0.6350 0.0780 0.1840]) 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.y(1:step:(simLen*10+1),15),'Color',[0 0.4470 0.7410]) 
    ylabel('HbO 5%CO_2','FontSize',fontS)
    subplot(rows,3,18)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),13),'Color',[0.4660 0.6740 0.1880]) %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),14),'Color',[0.6350 0.0780 0.1840]) 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.y(1:step:(simLen*10+1),15),'Color',[0 0.4470 0.7410]) 
    ylabel('HbO 10%CO_2','FontSize',fontS)
    figure(3)
    subplot(1,3,1)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),1),'r') %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),1),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),1),'b') 
    ylabel('N_NO','FontSize',fontS)
    subplot(1,3,2)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),2),'r') %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),2),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),2),'b') 
    ylabel('N_NPY','FontSize',fontS)
    subplot(1,3,3)
    plot(solC0.t(1:step:(simLen*10+1))+30,solC0.x(1:step:(simLen*10+1),3),'r') %NO NPY
    hold on
    plot(solC0.t(1:step:(simLen*10+1))+30,solC5.x(1:step:(simLen*10+1),3),'g') 
    plot(solC0.t(1:step:(simLen*10+1))+30,solC10.x(1:step:(simLen*10+1),3),'b') 
    ylabel('N_Pyr','FontSize',fontS)
end