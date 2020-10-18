% Function to simulate a linear transition from conditions in matrix M0 to
% conditions in matrix M1, over the period of time specified by 'interval.
% Can be used, for example, to simulate the gradual implementation of a lockdown
% over 10 days, if M0 gives pre-lockdown conditions, and M1 gives post-lockdown conditions 

function out = goveqs_scaleup(t, in, M0, M1, i, s, p, r, agg, sel, interval)

scale = min(max((t-interval(1))/(interval(2)-interval(1)),0),1);

M1.lin = M0.lin + scale*(M1.lin - M0.lin);
M1.lam = M0.lam + scale*(M1.lam - M0.lam);

out = goveqs_basis3(t, in, M1, i, s, p, r, agg, sel);