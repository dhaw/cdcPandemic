function f=phiPeaks
R0=2;
numseed=4;
lag=100:2.5:200;
llag=length(lag);
c=cell(llag,1);
%Feed in phi1&2?
load('forMAhpc.mat','C','Qeven')
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Qeven,0,0,R0);
maxPref=zeros(llag,1);
for i=1:llag
    lagi=lag(i);
    [~,~,~,tt,yy]=runPandemic2(gamma,n,nbar,na,NN,NNbar,NNrep,Kbar,K1,Cbar,Cbar,betaS,betaS,betaI,betaI,betaD,betaD,beta3,beta3,1,2,numseed,lagi,0);
    c{i}=[tt,yy];
    maxPref(i)=max(yy);
end
fs=12; lw=2;%30;%Font size
figure
cmap=parula(llag);
hold on
for i=1:llag
    ci=c{i};
    ti=ci(:,1); yi=ci(:,2); lti=length(ti);
    plot3(lag(i)*ones(lti,1),ti,yi,'-','linewidth',lw,'color',cmap(i,:))
end
hold off
xlabel('Lag')
ylabel('Time (days)')
zlabel('Prevalence')
axis ([lag(1),lag(end),1,500,0,.15])%max(maxPref)]);
%axis tight
title(strcat('R_0^{peak}=',num2str(R0)))
set(gca,'FontSize',fs);
grid on
grid minor
box on
