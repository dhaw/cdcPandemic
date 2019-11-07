function f=closureTimes(params)
close=7*(0:8);
lc=length(close);
c=cell(lc,1);
basePlot=pandemic1Dall(params,0,0,0,10^4);
tswitch=params(2);
for i=1:lc
    tclosei=tswitch+close(i);
    ploti=pandemic1Dall(params,0,0,0,tclosei);
    c{i}=ploti;%(ploti(:,1)>=tclosei,:);
end

figure
fs=12; lw=2;
%maxf=max([fall;ydata]);
hold on
cmap=flipud(parula(lc));
for i=1:lc
    vi=c{i};
    plot(vi(:,1),vi(:,2),'--','linewidth',lw,'color',cmap(i,:));
end
plot(basePlot(:,1),basePlot(:,2),'linewidth',lw,'color','k');
xlabel('Time (days)','FontSize',fs);
ylabel('Incidence','FontSize',fs);
set(gca,'FontSize',fs);
axis([0,360,0,max(basePlot(:,2))])
%legend([h1,h2],'Sim','Data','location','NW')
grid on
grid minor
box on
hold off