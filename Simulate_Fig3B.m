% v1: Version to iterate over different testing scenarios

clear all; Setup_parameters;
thresh1=1;
thresh2=220;%1500
xbeta=.5;
R0 = 2.5;
r.beta = r.beta*R0;

% --- Run the model -------------------------------------------------------
init = zeros(1,i.nx); seed = 10;
tmp = prm.N'; tmp2 = tmp(:)';
init(s.S) = tmp2;
init(i.IS1.urb.ad) = seed; init(i.S.urb.ad) = init(i.S.urb.ad) - seed;


% --- Baseline (no intervention)
tf = 2000;
M0 = make_model3(p, r, i, s, gps, prm);

geq = @(t,in) goveqs_basis3(t, in, M0, i, s, p, r, agg, sel);
[t,soln0] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx, 'RelTol', 1e-9, 'AbsTol', 1e-9));

% Find date of first death
vec = sum(soln0(:,i.aux.inc),2);%Mencodino %sum(soln1(:,i.aux.mort),2);
SIP0 = find(vec>1,1,'first');%Mendocino - SIP on day of first case %find(vec>1,1,'first');

% --- With lockdown over 10 days from date of first death
r1 = r; p1 = p;
r1.beta = r.beta*xbeta;
M1 = make_model3(p1, r1, i, s, gps, prm);

geq = @(t,in) goveqs_scaleup(t, in, M0, M1, i, s, p1, r, agg, sel, SIP0+[0 10]);
[t,soln1] = ode15s(geq, [0:1:tf], init, odeset('Nonnegative',1:i.nx, 'RelTol', 1e-9, 'AbsTol', 1e-9));

rq_list = linspace(1,14,25); %rq_list = 10;
pq_list = linspace(0,1,24);  %pq_list = 0.5;

for ir = 1:length(rq_list)
    fprintf('%0.5g ', ir);
    for ip = 1:length(pq_list)

        % --- Expansion of testing capacity
        t0_test  = 250;
        r2 = r1; p2 = p1;
        r2.q = 1/rq_list(ir);
        p2.q = pq_list(ip);
        M2 = make_model3(p2, r2, i, s, gps, prm);
        
        geq = @(t,in) goveqs_scaleup(t, in, M1, M2, i, s, p1, r, agg, sel, t0_test+[0 15]);
        [t,tmp] = ode15s(geq, [t0_test:1:tf], soln1(t0_test,:), odeset('Nonnegative',1:i.nx, 'RelTol', 1e-9, 'AbsTol', 1e-9));
        soln2 = soln1; soln2(t0_test+1:end,:) = tmp;
        
        % --- Relaxing of lockdown
        t0_relax = 400;
        r3 = r2; p3 = p2;
        r3.beta = r.beta;
        M3 = make_model3(p3, r3, i, s, gps, prm);
        
        geq = @(t,in) goveqs_scaleup(t, in, M2, M3, i, s, p1, r, agg, sel, t0_relax+[0 15]);
        [t,tmp] = ode15s(geq, [t0_relax:1:tf], soln2(t0_relax,:), odeset('Nonnegative',1:i.nx, 'RelTol', 1e-9, 'AbsTol', 1e-9));
        soln3 = soln2; soln3(t0_relax+1:end,:) = tmp;
        
        % Get the secondary peak
        t = 1:size(soln3,1); inds = find(t>t0_relax + 15);
        Hsol = squeeze(sum(soln3(:,s.H),2));
        Hmax(ir,ip) = max(Hsol(inds,:));
                
        allsol = cat(3,soln0, soln1, soln2, soln3);
        y = squeeze(sum(allsol(:,s.H,:),2));
        mat = diff(soln3,1); 
        num = mat(:,i.aux.hosp2);
        den = sum(mat(:,i.aux.hosp),2);
        
        vec = abs(diff(num./den)); vec(1:t0_relax) = 1;
        ind = find(vec<1e-5,1,'first');
        lvl(ir,ip) = num(ind)/den(ind);
    end
end
fprintf('\n');

figure; lw = 2; fs = 12;
h = contour(pq_list, rq_list, Hmax, [thresh1 thresh2], 'Color', 'b', 'linewidth', lw); 
xlabel('Proportion symptomatics diagnosed');
ylabel({'Average symptomatic duration','before diagnosis (days)'});
set(gca,'fontsize',fs);

grid on
box on

save scansurf_results;