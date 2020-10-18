clear all; 

ff=figure; fs = 14; lw = 2;

load illus_testing;

subplot(1,2,1);
y = squeeze(sum(allsol(:,s.H,:),2));
plot(y(:,3:5), 'linewidth', lw);
xlim([500 1500]);
line(xlim, 1416*[1 1], 'linestyle', ':', 'Color', 'k');
line(xlim, 916*[1 1], 'linestyle', ':', 'Color', 'k');
set(gca,'fontsize',fs);
xlabel('Days');
ylabel('Total number needing hospitalisation (all ages)');
title('Post-SIP, "resurgent" epidemic');
legend('No testing programme when SIP lifted','Modest testing programme','Aggressive testing programme');

load scansurf_results;
subplot(1,2,2); hold on;
h = contour(pq_list, rq_list, Hmax, [1 1500], 'Color', 'b', 'linewidth', lw); 
xlabel('Proportion symptomatics diagnosed');
ylabel({'Average symptomatic duration','before diagnosis (days)'});
set(gca,'fontsize',fs);

plot([0.25 0.75], [7, 4], '*', 'linewidth', 3, 'markersize', 12);
set(ff,'Position',[680   459   762   519]);