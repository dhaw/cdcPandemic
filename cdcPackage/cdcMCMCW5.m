function [xsto, outsto, history, accept_rate]=cdcMCMCW5(ydataNX)
% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

%MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);

n=10000;
sigma=1;
fixinds=[];
blockind=[];
displ=false;

%plim1=[3,6];
%plim2=[.3,.1];
plim1=[1.2,1.6];
plim2=[.25,.5];
x0=[randic(plim1),randic(plim2)];%,randic(plim3)];%,randic(plim4)];

F=@(params)fcn(params,ydataNX,plim1,plim2);%,plim3);%,plim4);
[xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
end

function f=unif(x,plim)
if (x-plim(1))*(x-plim(2))<0
    f=1/(plim(2)-plim(1));
else
    f=0;
end
end

function f=randic(plim)
f=plim(1)+rand*(plim(2)-plim(1));
end

function f=fcn(params,ydataNX,plim1,plim2)%,plim3)%,plim4)
f=-cdcLhoodsW5(params,ydataNX)+log(unif(params(1),plim1))+log(unif(params(2),plim2));%+log(unif(params(3),plim3));%+log(unif(params(4),plim4)));
end