function f=cdcDecision(xsto,ydata)%(pHN,int,thresh)%f and g from plotCumProb2
%Given input curve
%For each threshold
%Close schools or don't - run model
%
aggMonth=1;
monthShift=-1;
%}
qmean=1./[143.4400,364.7100,148.2000,64.6700];
NNbar=[19169690;62121035;184015269;39570590];
NN=sum(NNbar);
burn=200;
inc=1000;%Increment - every nth entry in MCMC chain
lx=length(xsto);
lout=floor((lx-burn)/inc);
c1=cell(lout,1);
c2=c1;
z=zeros(lout,2);
ind=1;
maxInc=0;
for i=burn+1:inc:lx 
    [~,out1,z1]=cdcPandemicSimulation(xsto(i,:),4:7,0,0,0,243);
    t1=out1(:,1); y1=out1(:,2:end); %z1=y1;%Output incidence by age
    %z1(t1<243,:)=[]; %t1(t1<243)=[];
    %z1=sum(z1,1).*qmean; z1=sum(z1);
    z1=sum(z1)/NN;
    [~,out2,z2]=cdcPandemicSimulation(xsto(i,:),4:7,0,0,0,360);
    t2=out2(:,1); y2=out2(:,2:end); %z2=y2;
    %z2(t2<243,:)=[]; %t2(t2<243)=[];
    %z2=sum(z2,1).*qmean; z2=sum(z2);
    z2=sum(z2)/NN;
    sumy1=sum(y1,2)/NN;
    sumy2=sum(y2,2)/NN;
    if aggMonth==1
       t1=1+ceil(t1/30+monthShift);
       sumy1=accumarray(t1,sumy1);
       t1=(t1(1):t1(end))';
       t2=1+ceil(t2/30+monthShift);
       sumy2=accumarray(t2,sumy2);
       t2=(t2(1):t2(end))';
    end
    maxInc=max(maxInc,max(sumy1));
    c1{ind}=[t1,sumy1];
    c2{ind}=[t2,sumy2];
    z(ind,:)=[z1,z2];
    ind=ind+1;
end
f=z;

fs=12; lw=.5;
cmap=lines(7);
figure
hold on
if aggMonth==1
    ydata=sum(ydata.*repmat(NNbar',size(ydata,1),1),2)/NN;
end
for i=1:lout
    c1i=c1{i};
    if aggMonth==1
        tswitch=9;
        h3=plot(1:length(ydata),ydata,'k-','linewidth',2);
        tvec=c1i(:,1);
        yvec=c1i(:,2);
        h1=plot(tvec(tvec<tswitch),yvec(tvec<tswitch),'-','color',cmap(1,:),'linewidth',lw);
        h2=plot(tvec(tvec>=tswitch-1),yvec(tvec>=tswitch-1),'-','color',cmap(2,:),'linewidth',lw);
    else
        tswitch=243;
        c2i=c2{i};
        tvec=c1i(:,1);
        yvec=c1i(:,2);
        h1=plot(tvec(tvec<tswitch),yvec(tvec<tswitch),'-','color',cmap(1,:),'linewidth',lw);
        h2=plot(tvec(tvec>=tswitch),yvec(tvec>=tswitch),'-','color',cmap(2,:),'linewidth',lw);
        %h1=plot(c1i(:,1),c1i(:,2),'-','color',cmap(1,:),'linewidth',lw);
        %h2=plot(c2i(:,1),c2i(:,2),'-','color',cmap(2,:),'linewidth',lw);
    end
end
%maxInc=max(c1i(:,2));
ylabel('Incidence')%('Hospitalisation')
%legend('0-4','5-19','2--64','65+')
set(gca,'fontsize',fs)
if aggMonth==1
    maxInc=max(maxInc,max(ydata));
    axis ([0,12,0,maxInc])
    xlabel('Time (months)')
    legend([h1,h2,h3],'Fit','Projection','Data','location','NW')
else
    axis ([0,360,0,maxInc])
    xlabel('Time (days)')
    legend([h1,h2],'Fit','Projection','location','NW')
end
grid on
grid minor
box on