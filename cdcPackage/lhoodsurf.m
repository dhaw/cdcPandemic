function f=lhoodsurf(NNbar,params,xdata,ydataNX)
p1=(.75:.001:.85);
p2=(0:.001:.04);
%
lp1=length(p1); lp2=length(p2);
prest=params(3:end);%mcmc=0
L=zeros(lp1,lp2);
for i=1:lp1
    p1i=p1(i);
    for j=1:lp2
        lhoodij=-sum(subLhoods(NNbar,[p1i,p2(j),prest],xdata,ydataNX));
        L(i,j)=lhoodij;
    end
end
f=L;
%}
figure
fs=10; lw=2;
imagesc(p1,p2,L')
hold on
%{
plot(.802,.0003,'o','markersize',5,'linewidth',1,'color','k','markerfacecolor','k')
text(.802+.005,.0003,'P_1','HorizontalAlignment','right','fontsize',12,'color','k')
plot(.78,.035,'o','markersize',5,'linewidth',1,'color','k','markerfacecolor','k')
text(.78+.005,.035,'P_2','HorizontalAlignment','left','fontsize',12,'color','k')
%}
set(gca,'YDir','normal')
set(gca,'fontsize',fs)
axis ([p1(1),p1(end),p2(1),p2(end)])
%caxis([-400,max(max(L))])
xlabel('x_2')
ylabel('\phi_2')
colorbar
grid on
grid minor
box on
hold off