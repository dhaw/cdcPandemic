function [xout,yout,NNbar]=subFSNextract(tab)
%tab=table of "extra" data - summer wave
%tab2=state=specific table of autumn wave
thisState='NM';%Match tab2

%if strcomp(tab.State(1),'State')==1
%    tab(1,:)=[];%If imported titles as a row
%en
%states=unique(tab.State);
%stateVec=strcmp(tab.State,thisState);
tab1=tab(tab.State==thisState,:);%(stateVec==1,:);
%tab1=tab(tab.State==thisState,:);
%Minnesota-like:
%
these={'State','agecat','MMWRWeek','weeklyrate','pop'};%cases/weeklyrate
tab1=tab1(:,ismember(tab1.Properties.VariableNames,these));
tab1=unstack(tab1,{'weeklyrate','pop'},'agecat');%cases/weeklyrate
yout=table2array(tab1(:,3:6));%3:7));%Ages
NNbar=table2array(tab1(:,7:end));%8:end));
NNbar=nanmean(NNbar,1);
NNmat=repmat(NNbar,size(yout,1),1);
yout=yout.*NNmat*10^(-5);
xout=tab1.MMWRWeek;
%}
%California-like:
%{
these={'CATCHMENT','AGECATEGORY','MMWRWEEK','WEEKLYRATE'};
tab1=tab1(:,ismember(tab1.Properties.VariableNames,these));
tab1=unstack(tab1,{'WEEKLYRATE'},'AGECATEGORY');%cases/weeklyrate
yout=table2array(tab1(:,3:7));
xout=tab1.MMWRWEEK;
NNbar=[];
%}
f=xout;

