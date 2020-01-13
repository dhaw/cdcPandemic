function f=profileLhoods(data,gamma)
param=(1.2:.01:1.7);
lp=length(param);
%Axis limits:
x0=factor(1);
xend=param(end);
%y0=-10;
%yend=120;
%
x=zeros(lp,1);
for i=1:lp
    x(i)=cdcLhoodsW5([param(i),gamma],data);
end
%For sum of lhoods:
%xsum=sum(x,1);
fs=12; lw=2;
figure
%colormap(parula(9))
hold on
%plot([x0,xend],[0,0],'k-','linewidth',lw)
%plot([1,1],[y0,yend],'k--','linewidth',lw)
plot(param,x,'-','linewidth',lw)
%plot(param,xsum,'k-','linewidth',lw)
hold off
xlabel('factor')
ylabel('-Log Likelihood')
axis tight%([x0,xend,y0,yend])
set(gca,'fontsize',fs')
grid on
grid minor
box on