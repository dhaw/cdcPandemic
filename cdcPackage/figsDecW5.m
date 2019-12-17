function [y1,y2,y3]=figsDecW5(ydata,ydataNX,xsto)%,y1,y2,y3)
NNbar=[1464566;3790730;10236474;3984200;2412129];
legString={'0-4','5-17','18-24','25-64','65+'};
nbar=length(NNbar);
%NNrep=repmat(NNbar',12,1);
%
    burn=1000;
    inc=100;
    lx=size(xsto,1);
    xdata=1:35;
    lxdata=length(xdata);
    nl=nbar*lxdata;
    Y=zeros(nl,floor((lx-burn)/inc));
    lz=size(Y,2);
    for i=1:lz%burn+1:inc:lx 
        ind=burn+(i-1)*inc+1;
        [yi,~,~]=cdcPandemicSimulationW5(xsto(ind,:),xdata,0,0,0,243);
        Y(:,i)=reshape(yi,nl,1);%.*NNrep,nl,1);
    end
    y123=prctile(Y,[2.3,50,97.5],2);
    y1=y123(:,1);
    y2=y123(:,2);
    y3=y123(:,3);
    %}

from=1;
to=16;
to2=35;
cmap=lines(7);
cmap2=cmap+.5*(1-cmap);
%y=ydata(1:35,:).*NNrep;
y=ydata;
maxy=max(max([y1;y2;y3]),max(max(y)));
fs=12; lw=2;
ticks=(1:8:35);
labs={'Apr','Jun','Aug','Oct','Dec'};
%
figure
hold on;
for j=1:nbar
    a1=lxdata*(j-1)+1; a2=lxdata*j;
    x1=y1(a1:a2);
    x2=y2(a1:a2);
    x3=y3(a1:a2);
    x=y(:,j);
    colj=cmap(j,:);
    colj2=cmap2(j,:);
    tvec=from:to;
    %
    plot(tvec,x1(tvec),'-','linewidth',1,'color',colj);
    plot(tvec,x2(tvec),'-','linewidth',1,'color',colj);
    tvec2= [tvec,fliplr(tvec)];
    inBetween=[x1(tvec);flipud(x3(tvec))];
    fill(tvec2,inBetween,colj,'facealpha',.2,'edgecolor',colj);
    plot(tvec,x2(tvec),'-','linewidth',lw,'color',colj);
    %}
    %Data:
    plot(tvec,x(tvec),'--','linewidth',lw,'color',colj);
    tvec=to:to2;
    %
    plot(tvec,x1(tvec),'-','linewidth',1,'color',colj2);
    plot(tvec,x2(tvec),'-','linewidth',1,'color',colj2);
    tvec2= [tvec,fliplr(tvec)];
    inBetween=[x1(tvec);flipud(x3(tvec))];
    fill(tvec2,inBetween,colj2,'facealpha',.2,'edgecolor',colj2);
    plot(tvec,x2(tvec),'-','linewidth',lw,'color',colj2);
    %}
    %Data:
    plot(tvec,x(tvec),'--','linewidth',lw,'color',colj);
    %}s
end
h1=plot([-1,-2],[-3,-4],'-','linewidth',lw,'color',cmap(1,:));
h2=plot([-1,-2],[-3,-4],'-','linewidth',lw,'color',cmap(2,:));
h3=plot([-1,-2],[-3,-4],'-','linewidth',lw,'color',cmap(3,:));
h4=plot([-1,-2],[-3,-4],'-','linewidth',lw,'color',cmap(4,:));
h5=plot([-1,-2],[-3,-4],'-','linewidth',lw,'color',cmap(5,:));
hold off
xlabel('Month')
ylabel('Incidence')
set(gca,'xtick',ticks,'xticklabel',labs)
set(gca,'fontsize',fs)
ax=gca;
ax.XAxis.MinorTickValues=tvec;
axis([xdata(1),xdata(end),0,maxy])%tvec(end)
legend([h1,h2,h3,h4,h5],legString,'location','NW')
grid on
grid minor
box on
