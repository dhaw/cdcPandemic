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

plim=[1.8,2.5;%min;max
    .25,.5]';
Cvec=[1.9200    0.4268    0.4600    0.4434    0.1665    1.7600    8.7522    2.4600    1.8630    1.2190    0.4500    0.4911    2.5900    0.8132    0.3412    4.5200    4.9913    5.6500    6.9995    4.3001    0.2300    0.2975    0.4300    0.6903    1.8930]';
plim=[[.9*Cvec,1.1*Cvec];plim]';

x0=randic(plim);

F=@(params)fcn(params,ydataNX,plim);
[xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, displ);
end

function f=unif(x,plim)
val=1./(plim(1,:)-plim(2,:));
in=(x-plim(1,:)).*(x-plim(2,:));
in(in>0)=0;
in(in<0)=1;
f=val.*in;
%{
if (x-plim(1))*(x-plim(2))<0
    f=1/(plim(2)-plim(1));
else
    f=0;
end
%}
end

function f=randic(plim)
f=plim(1,:)+rand(1,size(plim,2)).*(plim(2,:)-plim(1,:));
end

function f=fcn(params,ydataNX,plim)
f=-cdcLhoodsW5(params,ydataNX)+sum(log(unif(params,plim)));
%+log(unif(params(1),plim1))+log(unif(params(2),plim2));%+log(unif(params(3),plim3));%+log(unif(params(4),plim4)));
end