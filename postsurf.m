function f=postsurf(xsto)
burn=10000;
p1=(.8:.001:.804);
p2=(0:.00001:.0006);
xsto=xsto(burn+1:end,:);

figure;
fs=10; lw=2;
h=hist2d(xsto(:,1),xsto(:,2),100,'tile','probability');
imagesc(p1,p2,h)

axis ([p1(1),p1(end),p2(1),p2(end)])
caxis([0,max(max(h))])
xlabel('x_2')
ylabel('\phi_2')
colorbar
grid on
grid minor
box on
hold off