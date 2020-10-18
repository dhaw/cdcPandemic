function f=postsurf(xsto)
burn=10000;
p1=(.5:.001:.68);
p2=(0:.00001:.01);
xsto=xsto(burn+1:end,:);

figure;
fs=10; lw=2;
%h=hist2d(xsto(:,1),xsto(:,2),100,'tile','probability');
%imagesc(p1,p2,log10(h))
scatter(xsto(:,1),xsto(:,2),'o','markerfacecolor',.5*[1,1,1],'markeredgecolor',[1,1,1],'markerfacealpha',.1)
set(gca,'YDir','normal')

axis ([p1(1),p1(end),p2(1),p2(end)])
%caxis([0,max(max(h))])
xlabel('x_2')
ylabel('\phi_2')
%colorbar
grid on
grid minor
box on
hold off