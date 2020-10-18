function f=plotData(xdata,ydata)
nbar=5;
legString={'0-4','5-17','18-49','50-64','65+'};
figure
fs=10; lw=2;
thresh=35;
tend=length(xdata);%19
cmap=lines(nbar);
maxf=max(max(ydata));
hold on  
plotVec=xdata(1:tend);%1:size(ydata,1);%t1:t1+size(ydata,1)-1;
plot([thresh,thresh],[0,maxf],'k--','linewidth',lw)
h=zeros(1,nbar);
h(1)=plot(plotVec,ydata(1:tend,1),'-','linewidth',lw,'color',cmap(1,:));
h(2)=plot(plotVec,ydata(1:tend,2),'-','linewidth',lw,'color',cmap(2,:));
h(3)=plot(plotVec,ydata(1:tend,3),'-','linewidth',lw,'color',cmap(3,:));
h(4)=plot(plotVec,ydata(1:tend,4),'-','linewidth',lw,'color',cmap(4,:));
h(5)=plot(plotVec,ydata(1:tend,5),'-','linewidth',lw,'color',cmap(5,:));
set(gca,'FontSize',fs);
axis([17,52,0,maxf]);%([simCut+1,size(fall,1),0,maxf]);%tend
xlabel('Time (weeks)','FontSize',fs);
ylabel('Incidence','FontSize',fs);
legend(h,legString,'location','NW')
grid on
grid minor
box on
hold off