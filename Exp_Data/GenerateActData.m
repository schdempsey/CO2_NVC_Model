clear
load('rawAct.mat');
time=rawAct(:,1);
C0=[1 0 0];
C5=([0.6350 0.0780 0.1840]+[1 0 0])/2;
C10=[0.6350 0.0780 0.1840];
AZT=[0 0.85 0];

rawNormo15=rawAct(:,2:8);
meanNormo15=mean(rawNormo15,2);
stdNormo15=std(rawNormo15,[],2)./sqrt(length(rawNormo15(1,:)));

raw5C15=rawAct(:,9:15);
mean5C15=mean(raw5C15,2);
std5C15=std(raw5C15,[],2)./sqrt(length(raw5C15(1,:)));

raw10C15=rawAct(:,16:22);
mean10C15=mean(raw10C15,2);
std10C15=std(raw10C15,[],2)./sqrt(length(raw10C15(1,:)));

rawNormo15AZT=rawAct(:,23:29);
meanNormo15AZT=mean(rawNormo15AZT,2);
stdNormo15AZT=std(rawNormo15AZT,[],2)./sqrt(length(rawNormo15AZT(1,:)));

raw10C3=rawAct(:,30:34);
mean10C3=mean(raw10C3,2);
std10C3=std(raw10C3,[],2)./sqrt(length(raw10C3(1,:)));

ActData=cell(6,2);
ActData{1,1}=time;
ActData{2,1}=meanNormo15;
ActData{2,2}=stdNormo15;
ActData{3,1}=mean5C15;
ActData{3,2}=std5C15;
ActData{4,1}=mean10C15;
ActData{4,2}=std10C15;
ActData{5,1}=meanNormo15AZT;
ActData{5,2}=stdNormo15AZT;
ActData{6,1}=mean10C3;
ActData{6,2}=mean10C3;

