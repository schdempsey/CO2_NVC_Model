clear
load('rawExp4.mat');
time=rawExp4(:,1);
rawKO=rawExp4(:,2:8);
meanKO=mean(rawKO,2);
stdKO=std(rawKO,[],2)./sqrt(length(rawKO(1,:)));

rawCon=rawExp4(:,9:15);
meanCon=mean(rawCon,2);
stdCon=std(rawCon,[],2)./sqrt(length(rawCon(1,:)));

rawTxCon=rawExp4(:,16:23);
meanTxCon=mean(rawTxCon,2);
stdTxCon=std(rawTxCon,[],2)./sqrt(length(rawTxCon(1,:)));


figure(3)
errorbar(time,meanKO,stdKO,'r')
hold on
errorbar(time,meanCon,stdCon,'b')
errorbar(time,meanTxCon,stdTxCon,'k')
