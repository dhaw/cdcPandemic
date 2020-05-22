function [y1,y2,y3]=figsDec(ydata,ydataNX,xsto,y1,y2,y3)
NNbar=[19169690;62121035;184015269;39570590];
NNrep=repmat(NNbar',12,1);
age=1;
%%
%{
%fig1, fig2:
from=4;
to=12;
NNrep=repmat(NNbar',size(ydata,1),1);
y=sum(ydata.*NNrep,2);
maxy=max(y);
tvec=from:to;
fs=12; lw=2;
ticks=(4:2:12);
labs={'Apr','Jun','Aug','Oct','Dec'};
%
figure
hold on
%plot([jans(j),jans(j)],[0,mx],'k--','linewidth',lw)
plot(tvec,y(tvec),'k-','linewidth',lw);
hold off
xlabel('Month')
ylabel('Incidence')
set(gca,'xtick',ticks,'xticklabel',labs)
set(gca,'fontsize',fs)
ax=gca;
ax.XAxis.MinorTickValues=tvec;
axis([1,12,0,maxy])%tvec(end)
grid on
grid minor
box on
%}
%%
if age==0
    %{
    qmean=1./[143.4400,364.7100,148.2000,64.6700];
    burn=1000;
    inc=100;
    lx=size(xsto,1);
    %Z=zeros(floor((lx-burn)/inc),1);
    %lz=length(Z);
    Y=zeros(12,floor((lx-burn)/inc));
    lz=size(Y,2);
    for i=1:lz%burn+1:inc:lx 
        ind=burn+(i-1)*inc+1;
        [yi,~,~]=cdcPandemicSimulation(xsto(ind,:),1:12,0,0,0,243);
        Y(:,i)=sum(yi.*NNrep,2);%Z(i)=sum(z2.*qmean);
    end
    %{
    Z=[(1:lz)',Z];
    Z=sortrows(Z,2);
    ind1=Z(round(.025*lz),1);
    ind2=Z(round(.5*lz),1);
    ind3=Z(round(.975*lz),1);
    p1=xsto(burn+(ind1-1)*inc+1,:);
    p2=xsto(burn+(ind2-1)*inc+1,:);
    p3=xsto(burn+(ind3-1)*inc+1,:);
    [y1,~,~]=cdcPandemicSimulation(p1,1:12,0,0,0,243);
    [y2,~,~]=cdcPandemicSimulation(p2,1:12,0,0,0,243);
    [y3,~,~]=cdcPandemicSimulation(p3,1:12,0,0,0,243);
    %}
    y123=prctile(Y,[2.3,50,97.5],2);
    y1=y123(:,1);
    y2=y123(:,2);
    y3=y123(:,3);
    %%
    %fig3, fig4, fig5:
    %{
    y1=sum(y1.*NNrep,2);
    y2=sum(y2.*NNrep,2);
    y3=sum(y3.*NNrep,2);
    %}
    %}
    from=4;
    to=8;
    to2=12;
    col1=0*[1,1,1];
    col2=[1,0,0];
    y=sum(ydata(1:12,:).*NNrep,2);
    maxy=max([y;y1;y2;y3]);
    tvec=from:to;
    fs=12; lw=2;
    ticks=(4:2:12);
    labs={'Apr','Jun','Aug','Oct','Dec'};
    %
    figure
    plot(tvec,y1(tvec),'-','linewidth',1,'color',col1);
    hold on;
    plot(tvec,y2(tvec),'-','linewidth',1,'color',col1);
    tvec2= [tvec,fliplr(tvec)];
    inBetween=[y1(tvec);flipud(y3(tvec))];
    fill(tvec2,inBetween,col1,'facealpha',.2);
    %plot(tvec,y1(tvec),'b-','linewidth',lw);
    plot(tvec,y2(tvec),'k-','linewidth',lw);
    %plot(tvec,y3(tvec),'b-','linewidth',lw);
    %Data:
    plot(tvec,y(tvec),'k--','linewidth',lw);
    %
    tvec=to:to2;
    plot(tvec,y1(tvec),'-','linewidth',1,'color',col2);
    hold on;
    plot(tvec,y2(tvec),'-','linewidth',1,'color',col2);
    tvec2= [tvec,fliplr(tvec)];
    inBetween=[y1(tvec);flipud(y3(tvec))];
    fill(tvec2,inBetween,col2,'facealpha',.2,'edgecolor',col2);
    %plot(tvec,y1(tvec),'b-','linewidth',lw);
    plot(tvec,y2(tvec),'-','linewidth',lw,'color',col2);
    %plot(tvec,y3(tvec),'b-','linewidth',lw);
    %Data:
    plot(tvec,y(tvec),'k--','linewidth',lw);
    %}

    hold off
    xlabel('Month')
    ylabel('Incidence')
    set(gca,'xtick',ticks,'xticklabel',labs)
    set(gca,'fontsize',fs)
    ax=gca;
    ax.XAxis.MinorTickValues=tvec;
    axis([4,12,0,maxy])%tvec(end)
    grid on
    grid minor
    box on
else
    %{
    burn=1000;
    inc=100;
    lx=size(xsto,1);
    Y=zeros(48,floor((lx-burn)/inc));
    lz=size(Y,2);
    for i=1:lz%burn+1:inc:lx 
        ind=burn+(i-1)*inc+1;
        [yi,~,~]=cdcPandemicSimulation(xsto(ind,:),1:12,0,0,0,243);
        Y(:,i)=reshape(yi.*NNrep,48,1);%****NUMBER OF AGE GROUPS****
    end
    y123=prctile(Y,[2.3,50,97.5],2);
    y1=y123(:,1);
    y2=y123(:,2);
    y3=y123(:,3);
    %}
    %%
    %fig3, fig4, fig5:
    from=4;
    to=8;
    to2=12;
    cmap=lines(7);
    cmap2=cmap+.5*(1-cmap);
    y=ydata(1:12,:).*NNrep;
    maxy=max(max([y1;y2;y3]),max(max(y)));
    fs=12; lw=2;
    ticks=(4:2:12);
    labs={'Apr','Jun','Aug','Oct','Dec'};
    %
    figure
    hold on;
    for j=1:4
        a1=12*(j-1)+1; a2=12*j;
        x1=y1(a1:a2);
        x2=y2(a1:a2);
        x3=y3(a1:a2);
        x=y(:,j);
        colj=cmap(j,:);
        colj2=cmap2(j,:);
        tvec=from:to;
        %{
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
    hold off
    xlabel('Month')
    ylabel('Incidence')
    set(gca,'xtick',ticks,'xticklabel',labs)
    set(gca,'fontsize',fs)
    ax=gca;
    ax.XAxis.MinorTickValues=tvec;
    axis([4,12,0,maxy])%tvec(end)
    legend([h1,h2,h3,h4],'0-4','5-19','20-64','65+','location','NW')
    grid on
    grid minor
    box on
end
