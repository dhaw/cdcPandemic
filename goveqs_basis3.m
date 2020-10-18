function out = goveqs_basis3(t, in, M, i, s, p, r, agg, sel)

invec = in(1:i.nstates);
out   = zeros(length(in),1);

% --- Conctruct the linear terms
out(1:i.nstates) = M.lin*invec;

% Get new infections
lam = M.lam*invec;
newinfs = lam.*invec(s.S);

out(s.S) = out(s.S) - newinfs(:);
out(s.E) = out(s.E) + newinfs(:);

% Implement deaths
morts = M.mortvec.*invec;
out(1:i.nstates) = out(1:i.nstates) - morts;

% Get the auxiliaries - DH: don't worry too much about these for our
% current purposes, but will need them when looking at the data for daily
% hospitalisations. Can discuss then, how these work
out(i.aux.inc)   = agg.inc*(sel.inc.*M.lin)*invec;
out(i.aux.hosp)  = agg.hosp*(sel.hosp.*M.lin)*invec;
out(i.aux.mort)  = agg.mor*morts; 
out(i.aux.hosp2) = sum((sel.hosp2.*M.lin)*invec);