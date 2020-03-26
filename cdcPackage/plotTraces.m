function f=plotTraces(x1,x2,x3)
index=1;%size(x1,2)-4;
burn=0;
Y=[x1(:,index),x2(:,index),x3(:,index)];
a=0;%Lower bound in mcmc
b=10;%Upper bound in mcmc

figure
hold on
plot([burn,burn],[a,b],'k--','linewidth',2)
plot(Y)
axis([1,size(x1,1),a,b])
grid on
grid minor
box on