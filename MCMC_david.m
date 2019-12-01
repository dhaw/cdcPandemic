function [xsto, outsto, history, accept_rate]=MCMC_david(yvec)
% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

%MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);

%x0=[3,.58];
x0=[1.39,1/2.6];
n=10000;
sigma=1;
fixinds=[];
blockind=[];
displ=false;

plim1=[1.34,1.4];
plim2=[.33,.5];

F=@(params)fcn(params,yvec,plim1,plim2);
[xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
end

function f=unif(x,plim)
if (x-plim(1))*(x-plim(2))<0
    f=1/(plim(2)-plim(1));
else
    f=0;
end
end

function f=fcn(params,yvec,plim1,plim2)
f=-lhoods_david(params,yvec)+log(unif(params(1),plim1))+log(unif(params(2),plim2));
end