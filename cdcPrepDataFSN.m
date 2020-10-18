function f=cdcPrepDataFSN(tab)
%data=import data as table - default import options
%Produces ydata and ydataNX

ylab='Incidence';
state='New Mexico';

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
these={'MMWRYEAR','MMWRWEEK','AGECATEGORY','WEEKLYRATE'};
tab=tab(:,ismember(tab.Properties.VariableNames,these));
tab=unstack(tab,'WEEKLYRATE','AGECATEGORY');

f=tab;
%{

table(table.MMWRYEAR==2008,:)=[];
table(table.MMWRYEAR==2010,:)=[];
table(table.MMWRWEEK<=17,:)=[];
table.TOTAL=sum(table{:,3:end},2);

f=table;

tend=length(table.MMWRWEEK);
cmap=lines(7);
col1=cmap(1,:);
fs=12; lw=2;

figure
hold on
for i=1:size(table,2)-2
    plot(1:tend,table.TOTAL,'-','linewidth',lw,'color',cmap(i,:))
end
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel(ylab)
title(state)
axis tight
grid on
grid minor
box on
%}