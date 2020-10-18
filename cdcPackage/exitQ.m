function f=exitQ(M)%(pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta)
R0=2.5;
thresh1=916;
thresh2=1416;
d=1:.05:5;%7;
x=1./d;
lx=length(x);
p=0:.01:.75;
odds=p./(1-p);
lp=length(p);
%{
M=zeros(lx,lp);
for i=1:lx
    pr.qnew=x(i);
    Mi=zeros(1,lp);
    for j=1:lp
        pr.oddsnew=odds(j);
        [~,hosp]=seasonal1subAgeCovid19(pr,n,nbar,na,NN,NNbar,NNrep,Dout,beta);
        Mi(j)=max(hosp);
    end
   M(i,:)=Mi;
end
f=M;
%}
%
figure
fs=10; lw=2;
dhGray=.5*[1,1,1];
imagesc(p,d,M)
hold on
[C,h]=contour(p,d,M,[thresh1,thresh2],'LineColor','w','linewidth',lw,'ShowText','on');
plot(.2,4,'o','markersize',5,'linewidth',1,'color','w','markerfacecolor','w')
text(.2-.04,4,'P_1','HorizontalAlignment','right','fontsize',12,'color','w')
plot(.3,4.5,'o','markersize',5,'linewidth',1,'color','w','markerfacecolor','w')
text(.3+.02,4.5,'P_2','HorizontalAlignment','left','fontsize',12,'color','w')
clabel(C,h,'FontSize',12,'Color','w');
set(gca,'YDir','normal')
%{
plot([d(1),d(end)],[thresh1,thresh1],'--','linewidth',lw,'color',dhGray)
plot([d(1),d(end)],[thresh2,thresh2],'--','linewidth',lw,'color',dhGray)
plot(d,M,'k-','linewidth',lw);
%}
set(gca,'fontsize',fs)
axis tight%([p(1),p(end),d(1),d(end)])
xlabel('Asymptomatics: proportion identified and isolated')%('Time to quarantine')%('Proportion traced')
ylabel('Symptomatics: days from onset to isolation')%('Peak Hosp.')%('Time to quarantine')
title(strcat('R_0=',num2str(R0)))
colorbar
grid on
grid minor
box on
hold off
%}
        