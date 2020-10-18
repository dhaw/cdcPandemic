%function f=makeScenarios(vaxIn,NNbar,params,xdata,ydata)
function [F,Z]=makeScenarios(vaxIn,NNbar,params,xdata,ydata,xstoFull)
%function [Func,Zunc]=makeScenarios(vaxIn,NNbar,params,xdata,ydata,xsto)
numScen=8;%8;%24;%Change
na=5;
numDays=52*7;%720-243;%Change
%{
t1=(1:52)';
t2=7*ones(52,1);
day2week=repelem(t1,t2);
%}
%
mu=[393.0256  999.3054  702.6867  406.0680  177.1958];
%Epi:
R0=params(end-1);%1.4907; %params(end-1);%[1.45,1.7];%Change
muTimes=[1*ones(5,1),[.5*ones(4,1);1]];%Divide incidence

%Vax:
vaxon=[0,1];%[0,1];%[1e6,107,157];%Administration %Change
vaxparams=vaxIn;%Change (comment out next for 2009)
%{
vaxparams=ones(5,3);
vaxparams(:,1)=[.2,.5,.2,.09,.01]'/7;%Rate 1
vaxparams(:,2)=[0,0,.25,.35,.4]'/7;%Rate 2
vaxparams(:,3)=vaxIn(:,1,2);%Efficacy
%}

%Sch:
dayx=[1e6,61];
trigger=[2,.01];
%Antivirals:
tau=[0,1];
%%
%Fixed:
t0=243;
vaxSwitch=1.25e8/sum(NNbar);
%%
%Scenario numbers:
vax1=[2,4,7,8];%[2,4,7,8];%[2,4,5,9,10,14,21,22];*Change
vax2=0;%0;%[3,6,7,11,12,15,23,24];%Change
sch1=[5,7];%[5,7];%[4,6,9,11,16,18,21,23];%Change
sch2=sch1+1;%[5,7];%sch1+1;
av=[3,4,7,8];%[3,4,7,8];%[8:12,20:24];%Change
%%
scenariosAll=cell(numScen,1);%Cell array of "scenarios"
for i=1:numScen
    sceni=struct;
    %Epi:
    if i<=numScen/2
        sceni.R0=R0(1);
        sceni.muTimes=muTimes(:,1);
    else
        sceni.R0=R0(1);%Change - 1 for 2009, 2 for hyp
        sceni.muTimes=muTimes(:,1);%Change - 1 for 2009, 2 for hyp
    end
    %Vax::
    if ismember(i,vax1)
        sceni.vaxon=vaxon(2);
    elseif ismember(i,vax2)
        sceni.vaxon=vaxon(3);
    else
        sceni.vaxon=vaxon(1);
    end
    sceni.vaxparams=vaxparams;
    %Sch:
    if ismember(i,sch1)
        sceni.dayx=dayx(1);
        sceni.trigger=trigger(2);
        sceni.tswitch=243;
    elseif ismember(i,sch2)
        sceni.dayx=dayx(2);
        sceni.trigger=trigger(1);
        sceni.tswitch=287;
    else
        sceni.dayx=dayx(1);
        sceni.trigger=trigger(1);
        sceni.tswitch=243;
    end
    %Antivirals:
    if ismember(i,av)
        sceni.tau=tau(2);
    else
        sceni.tau=tau(1);
    end
    sceni.tv=181;%243;%181;%Change
    sceni.av=[0,0,.37*.29,.43*.2,.56*.17]';
    scenariosAll{i}=sceni;
end
f=scenariosAll;
%}
%%
%
%Single simulations
F=cell(1,numScen);
Z=zeros(numScen,na);
xsto=xstoFull(1001:100:end,:);
lx=size(xsto,1);
for i=1:numScen
    Fi=zeros(numDays,20,lx);
    for j=1:lx
        %[f,~,z2]=subPandemicSimulationVaxScenarios(NNbar,params,xdata,0,0,ydata,243,scenariosAll{i});
        [f,~,z2]=subPandemicSimulationVaxScenarios09(NNbar,xsto(j,:),xdata,0,0,ydata,243,scenariosAll{i});
        %f=accumarray(day2week,f,1);
        Fi(:,:,j)=f(1:numDays,:);
        Z(i,:)=z2;
    end
    F{i}=Fi;
end
%}
%%
%{
%Uncertainty
lx=size(xsto,1);
burn=10000;
int=100;
xsto=xsto(burn+1:int:end,:);
lf=size(xsto,1);%floor(lx-burn)/int;
Func=cell(1,numScen);
Zunc=zeros(numScen,lf);
%
r0=xsto(:,end-1);
r0med=median(r0);
r0=r0*1.45/r0med;
xsto(:,end-1)=r0;
parfor i=1:numScen/2
    Fi=zeros(90,lf);%Manually set number
    Zi=zeros(1,lf);
    sceni=scenariosAll{i};
    for j=1:lf
        [f,~,z2]=subPandemicSimulationVaxScenarios(NNbar,xsto(j,:),xdata,0,0,ydata,243,sceni);
        Fi(:,j)=sum(f,2);
        Zi(j)=sum(z2);
    end
    Func{i}=Fi;
    Zunc(i,:)=Zi;
end
r0=xsto(:,end-1);
r0med=median(r0);
r0=r0*1.7/r0med;
xsto(:,end-1)=r0;
parfor i=numScen/2+1:numScen
    Fi=zeros(90,lf);%Manually set number
    Zi=zeros(1,lf);
    for j=1:lf
        sceni=scenariosAll{i};
        [f,~,z2]=subPandemicSimulationVaxScenarios(NNbar,xsto(j,:),xdata,0,0,ydata,243,sceni);
        Fi(:,j)=sum(f,2);
        Zi(j)=sum(z2);
    end
    Func{i}=Fi;
    Zunc(i,:)=Zi;
end
%}