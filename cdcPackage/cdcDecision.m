function f=cdcDecision(xsto)%(pHN,int,thresh)%f and g from plotCumProb2
%Given input curve
%For each threshold
%Close schools or don't - run model
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
    [~,out1,~]=cdcPandemicSimulation(xsto(i,:),4:7,0,0,0,243);
    t1=out1(:,1); y1=out1(:,2:end); z1=y1;%Output incidence by age
    z1(t1<243,:)=[]; %t1(t1<243)=[];
    z1=sum(z1,1).*qmean; z1=sum(z1);
    [~,out2,~]=cdcPandemicSimulation(xsto(i,:),4:7,0,0,0,360);
    t2=out2(:,1); y2=out2(:,2:end); z2=y2;
    z2(t2<243,:)=[]; %t2(t2<243)=[];
    z2=sum(z2,1).*qmean; z2=sum(z2);
    sumy1=sum(y1,2)/NN;
    maxInc=max(maxInc,max(sumy1));
    c1{ind}=[t1,sumy1];
    c2{ind}=[t2,sum(y2,2)/NN];
    z(ind,:)=[z1,z2];
    ind=ind+1;
end
f=z;

fs=12; lw=.5;
cmap=lines(7);
figure
hold on
for i=1:lout
    c1i=c1{i};
    c2i=c2{i};
    plot(c1i(:,1),c1i(:,2),'-','color',cmap(1,:),'linewidth',lw);
    plot(c2i(:,1),c2i(:,2),'k--','linewidth',lw);
end
%maxInc=max(c1i(:,2));
xlabel('Time (days)')
ylabel('Incidence')%('Hospitalisation')
%legend('0-4','5-19','2--64','65+')
set(gca,'fontsize',fs)
axis ([0,360,0,maxInc])
grid on
grid minor
box on