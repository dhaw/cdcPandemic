function f=switchPeaks
R0=2;
numseed=4;
xwitch=100:5:300;
lxwitch=length(xwitch);
c=cell(lxwitch,1);
%Feed in phi1&2?
load('forMAhpc.mat','C','Qeven')
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Qeven,0,0,R0);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,CbarB,betaSb,betaIb,betaDb,beta3b,ages0]=prepFluAgeLocsB(C,Qeven,0,0,R0);
maxPrev=zeros(lxwitch,1);
for i=1:lxwitch
    xwitchi=xwitch(i);
    [~,~,~,tt,yy]=runPandemic2(gamma,n,nbar,na,NN,NNbar,NNrep,Kbar,K1,Cbar,CbarB,betaS,betaSb,betaI,betaIb,betaD,betaDb,beta3,beta3b,1,2,numseed,0,xwitchi);
    c{i}=[tt,yy];
    maxPrev(i)=max(yy);
end
fs=12; lw=2;%30;%Font size
figure
cmap=parula(lxwitch);
hold on
for i=1:lxwitch
    ci=c{i};
    ti=ci(:,1); yi=ci(:,2); lti=length(ti);
    plot3(xwitch(i)*ones(lti,1),ti,yi,'-','linewidth',lw,'color',cmap(i,:))
end
hold off
xlabel('Switch time')
ylabel('Time (days)')
zlabel('Prevalence')
axis ([xwitch(1),xwitch(end),1,500,0,.15])%max(maxPrev)]);
%axis tight
title(strcat('R_0^{term}=',num2str(R0)))
set(gca,'FontSize',fs);
grid on
grid minor
box on
