% DH: This code simulates both the lockdown and the release. To plot the
% panels in Fig.1, keep ylim as [0 400] as in line 61. To plot Figure 2,
% change that to ylim([0 tf])

%clear all; 

Setup_parameters;

thresh1=115;%1416
thresh2=220;%916
xlimupper=1000;
ylimupper=2200;
xbeta=.75;

R0 = 2.5;
r.beta = r.beta*R0;

% --- Run the model -------------------------------------------------------
init = zeros(1,i.nx); seed = 10;
tmp = prm.N'; tmp2 = tmp(:)';
init(s.S) = tmp2;
init(i.IS1.urb.ad) = seed; init(i.S.urb.ad) = init(i.S.urb.ad) - seed;


tf = 1000;

% --- 1. Baseline (no intervention)
M0 = make_model3(p, r, i, s, gps, prm);

geq = @(t,in) goveqs_basis3(t, in, M0, i, s, p, r, agg, sel);
[t,soln1] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));

% Find date of first death
vec = sum(soln1(:,i.aux.inc),2);%Mencodino %sum(soln1(:,i.aux.mort),2);
t_ldown = find(vec>1,1,'first');%Mendocino - SIP on day of first case %find(vec>1,1,'first');
t_relax = tf/2;


% --- 2. With lockdown over 10 days from date of first death
r1 = r; p1 = p;
r1.beta = r.beta*xbeta;
M1 = make_model3(p1, r1, i, s, gps, prm);

geq = @(t,in) goveqs_scaleup(t, in, M0, M1, i, s, p1, r, agg, sel, t_ldown+[0 10]);
[t,tmp1] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));

% Then lifting the lockdown
geq = @(t,in) goveqs_basis3(t, in, M0, i, s, p, r, agg, sel);
[t,tmp2] = ode15s(geq, [t_relax:1:tf], tmp1(t_relax,:), odeset('Nonnegative',1:i.nx));

soln2 = tmp1; soln2(t_relax+1:end,:) = tmp2;



% --- Plot the results ----------------------------------------------------

figure; fs = 12; lw = 1.5;
hold on;

allsol = cat(3, soln1, soln2);

y = allsol(:,intersect(s.urb,s.H),:) + allsol(:,intersect(s.rur,s.H),:);

p1 = plot(y(:,:,1),'-', 'linewidth', lw);
set(gca,'ColorOrderIndex',1);
plot(y(:,:,2),'--', 'linewidth', lw);
% set(gca,'ColorOrderIndex',1);
% plot(y3,':', 'linewidth', lw);
 xlim([0 xlimupper]);
 ylim([0,ylimupper]);

line(xlim, thresh1*[1 1], 'linestyle', ':', 'Color', 'k');
line(xlim, thresh2*[1 1], 'linestyle', ':', 'Color', 'k');

legend(p1, '<15yo', '16 - 64yo', '>65yo');

set(gca,'fontsize',fs);
xlabel('Time (days)')
ylabel('Hospitalisations (25% reduction)')
title(strcat('R_0=',num2str(R0)))
box on