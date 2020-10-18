function f=makeFigMLE(f1,f2,f3,xdata,ydata)
Y=[sum(f1,2),sum(f2,2),sum(f3,2)];
Ydata=sum(ydata,2);
maxy=max(max(max(Y)),max(Ydata));
fs=10; lw=2;
cmap=lines(7);
one3=ones(1,3);
col1=cmap(1,:);%.8*one3;
col2=cmap(2,:);%.4*one3;
col3=cmap(3,:);%0*one3;
figure
hold on
plot(xdata,Y(:,1),'linewidth',lw,'color',col1)
plot(xdata,Y(:,2),'linewidth',lw,'color',col2)
plot(xdata,Y(:,3),'linewidth',lw,'color',col3)
plot(xdata,Ydata,'k--','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel('Incidence')
axis([xdata(1),xdata(end),0,maxy])
legend('C_{22}','\beta_2','C (full)','Data','location','NW')
grid on
grid minor
box on