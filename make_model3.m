% DH: This function creates full set of matrices for a given set of
% parameters 'r' (for rates) and 'p' (for proportions)

function M = make_model(p, r, i, s, gps, prm)

m = zeros(i.nstates);

% DH: All diagonal elements are positive, and off-diagonal are negative.
% Any element in row i and column j (i not equal to j) designates the
% per-capita rate from compartment j to compartment i

for ic = 1:length(gps.geo)
    geo = gps.geo{ic};
    
    for ia = 1:length(gps.age)
        age = gps.age{ia};
        
        S   = i.S.(geo).(age);       % Uninfected
        E   = i.E.(geo).(age);       % Exposed
        IA  = i.IA.(geo).(age);      % Asymptomatic
        IP  = i.IP.(geo).(age);      % Pre-symptomatic, mild
        IN1 = i.IN1.(geo).(age);     % Symptomatic, to be missed by testing strategy
        IN2 = i.IN2.(geo).(age);     
        IS1 = i.IS1.(geo).(age);     % Symptomatic, to be identified by testing strategy
        IS2 = i.IS2.(geo).(age);     
        Q1  = i.Q1.(geo).(age);      % Quarantine
        Q2  = i.Q2.(geo).(age);
        H   = i.H.(geo).(age);       % Hospitalization
        R   = i.R.(geo).(age);       % recovered
        
        % Incubation
        source  = E;
        destins =           [IA,           IP];
        rates   = [(1-p.sympto),     p.sympto]*r.incub;
        m(destins, source) = m(destins, source) + rates';
        
        % DH: See note above. Lines 32-35, for example, are setting up a
        % transition from E to the symptomatic and pre-symptomatic
        % compartments
        
        
        % Pre-symptomatic to symptomatic
        source  = IP;
        destins =                    [IN1,                 IN2,                  IS1,              IS2];
        rates   = [(1-p.q)*(1-p.hosp(ia)),  (1-p.q)*p.hosp(ia),   p.q*(1-p.hosp(ia)),   p.q*p.hosp(ia)]*r.eta;
        m(destins, source) = m(destins, source) + rates';
        
        % Quarantine
        source  = IS1;
        destin  = Q1;
        rate    = r.q;
        m(destin, source) = m(destin, source) + rate;

        source  = IS2;
        destin  = Q2;
        rate    = r.q;
        m(destin, source) = m(destin, source) + rate;
        
        % Hospitalization
        sources = [IS2, IN2, Q2];
        destin  = H; 
        rates   = r.hosp;
        m(destin, sources) = m(destin, sources) + rates;
        
        % --- Recovering from disease
        sources = [IA IS1 IN1 Q1];
        destin  = R;
        rate    = r.gamma;
        m(destin, sources) = m(destin, sources) + rate;
        
        sources = [Q2 H];
        destin  = R;
        rate    = r.gamma_h(ia);
        m(destin, sources) = m(destin, sources) + rate;
    end
    
end

M.lin = sparse(m - diag(sum(m,1)));
% DH: M.lin is now the matrix containing the linear (non-force-of-infection) transitions 


% --- Set up the lambda matrices

% The motif to be replicated for each infectious class
motif = repmat(prm.mixmat, 1, length(gps.geo));
% Adjust motif by population numbers
tmp = prm.N'; tmp2 = tmp(:)'; 
motif = motif./repmat(tmp2,length(gps.age),1);

% Now implement the motif to construct the set of rows to be used in the
% FOI, for each geography
block = zeros(length(gps.age), i.nstates);
block(:,s.IA)  = p.c*motif;
block(:,s.IP)  = p.c*motif;
block(:,s.IN1) = motif;
block(:,s.IN2) = motif;
block(:,s.IS1) = motif;
block(:,s.IS2) = motif;
% Discount infection in rural settings
block(:,s.rur) = block(:,s.rur)*p.rurbeta;

% Implement inter-patch connectivity
lg = length(gps.geo); la = length(gps.age);
ipmat = zeros(la*lg);
for ir = 1:lg
    rows = (ir-1)*la + [1:la];
    for ic = 1:lg
        cols = (ic-1)*la + [1:la];
        ipmat(rows,cols) = prm.connmat(ir,ic);
    end
end
M.lam = r.beta*repmat(block,lg,1).*repmat(ipmat,1,i.nstates/(lg*la));

% DH: M.lam is now a matrix that will be multiplied by the state vector to
% give the separate forces-of-infection for children, adults and elderly,
% in both urban and rural settings


% --- Get the mortality rates
m = zeros(1,i.nstates);
m(intersect(s.H, s.ch)) = r.mu(1);                                     % Disease induced mortality (children)
m(intersect(s.H, s.ad)) = r.mu(2);                                     % Disease induced mortality (young adult)
m(intersect(s.H, s.el)) = r.mu(3);                                    % Disease induced mortality (old)

M.mortvec = m';
