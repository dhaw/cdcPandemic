function f=profileLhoods(data)
factor=.005:.001:1;
lx=4;
%Axis limits:
x0=factor(1);
xend=factor(end);
y0=-10;
yend=120;
%
lf=length(factor);
dim=4*lx;%length(xdata)*4
x=zeros(dim,lf);
for i=1:lf
    x(:,i)=lhoods(data,factor(i));
end
%For sum of lhoods:
xsum=sum(x,1);
fs=12; lw=2;
figure
colormap(parula(9))
hold on
plot([x0,xend],[0,0],'k-','linewidth',lw)
plot([1,1],[y0,yend],'k--','linewidth',lw)
plot(factor,x,'-','linewidth',1)
plot(factor,xsum,'k-','linewidth',lw)
hold off
xlabel('factor')
ylabel('Log Likelihood')
axis([x0,xend,y0,yend])
set(gca,'fontsize',fs')
grid on
grid minor
box on