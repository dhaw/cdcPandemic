function [y1,y2,y3]=figsDecW5data(ydata,mmean,msd)
NNbar=[1464566;3790730;10236474;3984200;2412129];
legString={'0-4','5-17','18-24','25-64','65+'};
nbar=length(NNbar);
lt=size(ydata,1);
NNrep=repmat(NNbar',lt,1);
mrep=repmat(mmean,lt,1);
msdrep=repmat(msd,lt,1);

x2=ydata.*mrep;
x1=ydata.*(mrep-msdrep);
x3=ydata.*(mrep+msdrep);
x4=sum(x2,2);

from=1;
to=16;
to2=35;
tvec=from:to2;
cmap=lines(7);
%cmap2=cmap+.5*(1-cmap);
%y=ydata(1:35,:).*NNrep;
maxy=max(max(max(x3)),max(x4));
fs=12; lw=2;
ticks=tvec(1:8:35);
labs={'Apr','Jun','Aug','Oct','Dec'};
%
figure
hold on;
for j=1:nbar
    colj=cmap(j,:);
    plot(tvec,x1(tvec,j),'--','linewidth',1,'color',colj);
    plot(tvec,x2(tvec,j),'-','linewidth',1,'color',colj);
    tvec2= [tvec,fliplr(tvec)];
    inBetween=[x1(tvec,j);flipud(x3(tvec,j))];
    fill(tvec2,inBetween,colj,'facealpha',.2,'edgecolor',colj);
    %plot(tvec,x2(tvec),'-','linewidth',lw,'color',colj);
end
plot(tvec,x4(tvec),'k-','linewidth',lw)
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
axis([1,length(tvec),0,maxy])%tvec(end)
legend([h1,h2,h3,h4,h5],legString,'location','NW')
grid on
grid minor
box on
