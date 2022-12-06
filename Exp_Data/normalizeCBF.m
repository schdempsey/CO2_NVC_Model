clc
clear
Name={'Exp1','Exp2','Exp3'};
for j=1:3
Expstruct=load(strcat('raw',Name{j},'.mat'));
Exp = Expstruct.(sprintf(Name{j}));
[rows,cols]=size(Exp);
normo=Exp{1,1};
baseline=mean(normo(1:4));
for i=1:(cols-1)
    for k=1:3
        Mat=Exp{k,i};
        Exp{k,i}=Mat./baseline;
    end
end
Data=Exp;
save(strcat(Name{j},'.mat'),'Data')
end