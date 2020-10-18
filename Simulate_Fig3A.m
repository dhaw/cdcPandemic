% v2: Slightly tidied up version, for use in writeup 

clear all; 

Setup_parameters;

thresh1=115;%1416
thresh2=220;%916
xlimupper=1400;
xbeta=.5;

R0 = 2.5; 
r.beta = r.beta*R0;

% --- Run the model -------------------------------------------------------
init = zeros(1,i.nx); seed = 10;
tmp = prm.N'; tmp2 = tmp(:)';
init(s.S) = tmp2;
init(i.IS1.urb.ad) = seed; init(i.S.urb.ad) = init(i.S.urb.ad) - seed;


% --- 0. Baseline (no intervention)
tf = 2000;

M0 = make_model3(p, r, i, s, gps, prm);

geq = @(t,in) goveqs_basis3(t, in, M0, i, s, p, r, agg, sel);
[t,soln0] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx));

% Find date of first death
vec = sum(soln0(:,i.aux.inc),2);%Mencodino %sum(soln1(:,i.aux.mort),2);
SIP0 = find(vec>1,1,'first');%Mendocino - SIP on day of first case %find(vec>1,1,'first');

t_ldown = SIP0;
t_lrelx = 500;


% --- 1. With lockdown over 10 days from date of first death
r1 = r; p1 = p;
r1.beta = r.beta*xbeta;
M1 = make_model3(p1, r1, i, s, gps, prm);

geq = @(t,in) goveqs_scaleup(t, in, M0, M1, i, s, p1, r, agg, sel, t_ldown+[0 10]);
[t,tmp] = ode15s(geq, [t_ldown:1:tf], soln0(t_ldown,:), odeset('Nonnegative',1:i.nx));
soln1 = soln0; soln1(t_ldown+1:end,:) = tmp;


% --- 2. With release of lockdown at time t_lrelx, no testing
r2 = r1; p2 = p1;
r2.beta = r.beta;
M2 = make_model3(p2, r2, i, s, gps, prm);

geq = @(t,in) goveqs_scaleup(t, in, M1, M2, i, s, p2, r2, agg, sel, t_lrelx+[0 10]);
[t,tmp] = ode15s(geq, [t_lrelx:1:tf], soln1(t_lrelx,:), odeset('Nonnegative',1:i.nx));
soln2 = soln1; soln2(t_lrelx+1:end,:) = tmp;


% --- 3. With release of lockdown at time t_lrelx, modest testing
r3 = r2; p3 = p2;
r3.q = 1/7; p3.q = 0.25;
M3 = make_model3(p3, r3, i, s, gps, prm);

geq = @(t,in) goveqs_scaleup(t, in, M1, M3, i, s, p3, r3, agg, sel, t_lrelx+[0 10]);
[t,tmp] = ode15s(geq, [t_lrelx:1:tf], soln1(t_lrelx,:), odeset('Nonnegative',1:i.nx));
soln3 = soln1; soln3(t_lrelx+1:end,:) = tmp;


% --- 4. With release of lockdown at time t_lrelx, aggressive testing
r4 = r2; p4 = p2;
r4.q = 1/4; p4.q = 0.75;
M4 = make_model3(p4, r4, i, s, gps, prm);

geq = @(t,in) goveqs_scaleup(t, in, M1, M4, i, s, p4, r4, agg, sel, t_lrelx+[0 10]);
[t,tmp] = ode15s(geq, [t_lrelx:1:tf], soln1(t_lrelx,:), odeset('Nonnegative',1:i.nx));
soln4 = soln1; soln4(t_lrelx+1:end,:) = tmp;

allsol = cat(3,soln0,soln1,soln2,soln3,soln4);

figure; fs = 12; lw = 2;

y = squeeze(sum(allsol(:,s.H,:),2));
plot(y(:,3:5), 'linewidth', lw);
xlim([500 xlimupper]);
line(xlim, thresh1*[1 1], 'linestyle', ':', 'Color', 'k');
line(xlim, thresh2*[1 1], 'linestyle', ':', 'Color', 'k');
set(gca,'fontsize',fs);
xlabel('Days');
ylabel('Total number needing hospitalisation (all ages)');
title('Post-SIP, "resurgent" epidemic');
legend('No testing programme when SIP lifted','Modest testing programme','Aggressive testing programme');

box on

save illus_testing;
