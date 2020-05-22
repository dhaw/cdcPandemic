function [xdata2,ydata2nx,ydata2]=subFSNextractState(tab2,NNbar,xdata,ydata)
%tab2=import data as table - default import options
%Rest - stare-specific output from subFSNextract
%Produces ydata and ydataNX for whole


ylab='Incidence';
state='New York';
legString={'0-4','5-17','18-49','50-64','65+'};

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
these={'MMWRYEAR','MMWRWEEK','AGECATEGORY','WEEKLYRATE'};
tab2=tab2(:,ismember(tab2.Properties.VariableNames,these));
tab2=unstack(tab2,'WEEKLYRATE','AGECATEGORY');
toHere=xdata(end);

tab2(tab2.MMWRYEAR==2008,:)=[];
tab2(tab2.MMWRYEAR==2010,:)=[];
tab2(tab2.MMWRWEEK<=toHere,:)=[];
%table.TOTAL=sum(table{:,3:end},2);
ydata2=table2array(tab2(:,3:end));

%Check per state
%Minnesota:%Incorrect download?
%ydata2=[ydata2(:,[1,3,2,4]),sum(ydata2(:,5:7),2)];%Incorrect download?
%California:
%Leave it!

NNmat=repmat(NNbar,size(ydata2,1),1);
ydata2=ydata2.*NNmat*10^(-5);

xdata2=tab2.MMWRWEEK;
xdata2(xdata2<=toHere)=[];

xdata2=[xdata;xdata2];
ydata2nx=[ydata;ydata2];
ydata2=ydata2nx.*repmat(mu,size(ydata2nx,1),1);
ydata2(isnan(ydata2)==1)=0;
%
cmap=lines(7);
col1=cmap(1,:);
fs=12; lw=2;
figure
hold on
plot(xdata2,ydata2,'-','linewidth',lw)%,'color',cmap(i,:))
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel(ylab)
title(state)
legend(legString,'location','NW');
axis tight
grid on
grid minor
box on
%}