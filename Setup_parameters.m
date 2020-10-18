%clear all; 

gps.geo = {'urb','rur'};                            
gps.age = {'ch','ad','el'};

states1 = {'S','E','IA','IP','IN1','IN2','IS1','IS2','Q1','Q2','R','H'};

[i, s, d, lim] = get_addresses({states1, gps.geo, gps.age}, [], [], [], 0);
d = char(d);

% DH: Note that we've now created the following:
% -- i is a Matlab structure containing the addresses of each compartment.
% E.g. type i.IA.urb.ad to find that the compartment number for
% asymptomatic infections who are urban and adults is 14
% -- s is a Matlab compartment containing all the addresses of given groups
% of compartments, e.g. type s.IA to find the addresses of all asymptomatic
% infections. And type intersect(s.IA, s.ch) to find all the child
% asymptomatic infections
% -- d is a list of each compartment


% --- Include the auxiliaries ---------------------------------------------

% DH: These 'auxiliaries' are used to count cumulative numbers over time,
% e.g. incidence. They're all quantities involving integration, so we set
% up 'auxiliary matrices' that we use to perform these integrations.
% However note that we don't really use these in the calculations (we only
% show prevalence in H rather than cumulative hospitalisations), so you can
% skip these for now

numcoms = length(gps.geo)*length(gps.age);

names = {  'inc',     'hosp',     'hosp2',   'mort',  'inc2'};
lgths = [numcoms,    numcoms,           1,  numcoms,       6];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;

s.infectious = [s.IA, s.IP, s.IN1, s.IN2, s.IS1, s.IS2];                         % Infectious
s.prevalent  = [s.IN1, s.IN2, s.IS1, s.IS2, s.H, s.Q1, s.Q2];   
 

% --- Counting symptomatic incidence, by age and city
mat = zeros(i.nstates);
mat([s.IN1, s.IN2, s.IS1,s.IS2],:) = 1;
sel.inc = sparse(mat - diag(diag(mat)));

mat = zeros(numcoms,i.nstates); row = 1;
for ic = 1:length(gps.geo)
    for ia = 1:length(gps.age)
        mat(row, [i.IN1.(gps.geo{ic}).(gps.age{ia}), i.IN2.(gps.geo{ic}).(gps.age{ia}), i.IS1.(gps.geo{ic}).(gps.age{ia}), i.IS2.(gps.geo{ic}).(gps.age{ia})]) = 1;
        row = row+1;
    end    
end
agg.inc = sparse(mat);  


% --- Counting deaths
mat = zeros(numcoms,i.nstates); row = 1;
for ic = 1:length(gps.geo)
    for ia = 1:length(gps.age)
        mat(row, [intersect(intersect(s.H,s.(gps.age{ia})),s.(gps.geo{ic}))]) = 1;
        row = row+1;
    end    
end
agg.mor = sparse(mat);  


% --- Counting hospitalisations
mat = zeros(i.nstates);
mat(s.H,:) = 1;
sel.hosp = sparse(mat - diag(diag(mat)));

mat = zeros(numcoms,i.nstates); row = 1;
for ic = 1:length(gps.geo)
    for ia = 1:length(gps.age)
        mat(row, i.H.(gps.geo{ic}).(gps.age{ia})) = 1;
        row = row+1;
    end    
end
agg.hosp = sparse(mat);    


% --- Counting hospitalisations specifically coming from Q
mat = zeros(i.nstates);
mat(s.H,s.Q2) = 1;
sel.hosp2 = sparse(mat - diag(diag(mat)));




% --- Set up the parameters -----------------------------------------------

% DH: all rates for the model are contained in 'r', all parameters are contained in 'p'

R0 = 1;                                                                    % Basic reproduction number
incub_period  = 5;                                                         % Average incubation period (days)
presym_period = 1;                                                         % Average pre-symptomatic period (days)
infec_period  = 5;                                                         % Average infectious period (days)

prop_hosp    = [0.04 0.278 0.599];                                         % Age-specific proportions needing hospitalisation (severe disease)
cfr          = [0.01, 0.3, 6.4]/100;                                       % Age-specific case fatality rates
mort_on_hosp = cfr./prop_hosp;                                             % Amongst those being hospitalised, what proportion would die
dur_in_hosp  = 10;                                                         % Average days in hospital


% --- Set up parameter values
r.incub   = 1/incub_period;                                                % Rate of progression to infectiousness
p.sympto  = 2/3;                                                           % Proportion developing symptomatic infection
p.c       = 2/3;                                                           % Relative infectiousness of asymptomatic vs symptomatic infection
r.eta     = 1/presym_period;                                               % Rate from pre-symptomatic to symptomatic infection
r.gamma   = 1/infec_period;
p.hosp    = prop_hosp;
r.hosp    = 1/5;                                                           % Amongst those with severe disease, average rate of progression to needing hospitalisation
r.mu      = mort_on_hosp*(1/dur_in_hosp);                                  % Rate of mortality amongst those needing hospitalisation
r.gamma_h = (1-mort_on_hosp)*(1/dur_in_hosp);                              % Rate of recovery amongst those needing hospitalisation


% Interventions
r.q       = 0;
p.q       = 0;                                                             % fraction of symptomatic that to be quarantined

load mixmat; prm.mixmat = mixmat;                                          

% --- Demographic parameters
popn        = [15702,52630,19656];%Mendocino %[80816, 321096, 98030];                                      
purb        = 0.7;% 0.326;%Mendocino %0.7;
prm.N       = [popn*purb; popn*(1-purb)];
p.crossg    = 0.05;                                                        % Infectiousness between I in one setting and S in another
prm.connmat = [1 p.crossg; p.crossg 1];                                    
p.rurbeta   = 0.75;% 0.85;%Mendocino 0.75;

% --- Find beta so that we get given value of R0 in Urban -----------------
r.beta = 1;
M      = make_model3(p, r, i, s, gps, prm);
Mlin   = M.lin;
Mlam   = M.lam;

% Pick out indices for the setting where we're interested in calibrating R0
I_inds = intersect([s.E, s.IA, s.IP, s.IN1, s.IN2, s.IS1, s.IS2],s.urb);

% Adjust for the number susceptible
tmp = prm.N'; tmp2 = tmp(:);
Mlam = Mlam.*repmat(tmp2,1,i.nstates);

% Construct F matrix
mat = zeros(i.nstates);
mat(s.E,:) = Mlam;
F = mat(I_inds, I_inds);

% Construct V matrix
m = -(Mlin - diag(M.mortvec));
V = m(I_inds, I_inds);

% Get R0
maxeig = max(eigs(F*inv(V)));
r.beta = 1/maxeig;                                                         % Temporary value for beta, to be replaced when R0 is specified
