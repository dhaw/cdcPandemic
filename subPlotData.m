function f=subPlotData(yx)
thresh=10;
figure
hold on
plot([0,36],[thresh,thresh],'k--')
plot(yx,'linewidth',1.5)
plot(sum(yx,2),'k-','linewidth',2)
axis tight
grid on
grid minor
box on
hold off