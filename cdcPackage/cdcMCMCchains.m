function f=cdcMCMCchains(xsto4)%,xsto4b)
%names={'m','x_{22}','R_0','\gamma'};
names={'R_0','\gamma'};
burn=1000;
figure;
 for i=1:2
     subplot(2,1,i);
     %plot(1001:length(xsto4),xsto4(1001:end,i));
     plot(burn+1:length(xsto4),[xsto4(burn+1:end,i)]);%,xsto4b(burn+1:end,i)]);
     axis tight;
     title(names{i});
 end
 %{
 subplot(2,2,2);
 plot(1001:length(xsto4),xsto4(1001:end,2));
 title('x_{22}'); axis tight; subplot(2,2,3);
 plot(1001:length(xsto4),xsto4(1001:end,3));
 title('R_0');
 axis tight; subplot(2,2,4);
 plot(1001:length(xsto4),xsto4(1001:end,4));
 title('\gamma'); axis tight;
 %}