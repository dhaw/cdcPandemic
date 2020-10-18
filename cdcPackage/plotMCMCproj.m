function f=plotMCMCproj(xsto)
fs=10; lw=2; ms=5;
figure
%scatter(xsto(:,1),xsto(:,2),'o','markeredgecolor','k','markerfacecolor','k','markerfacealpha',.2)
hist2d(xsto(:,1),xsto(:,2),50)
xlabel('x_{22}')
ylabel('\phi_2')
set(gca,'fontsize',fs)
axis tight
grid on
grid minor
box on
