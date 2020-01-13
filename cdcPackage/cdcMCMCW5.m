function [xsto, outsto, history, accept_rate]=cdcMCMCW5(ydataNX,thetac)
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

plim=[1,3;%min;max
    .25,.5];
Cvec=[1.9472    0.4316    0.4571    0.4432    0.1658    1.9350    8.7949    2.4030    1.8097    1.2314    0.4651    0.4938    2.3491    0.8098    0.3408    4.5977    5.0972    5.4579  6.9881    4.2751    0.2365    0.2976    0.4186    0.7152    1.8597]';
plim=[[.9*Cvec,1.1*Cvec];plim]';

%x0=randic(plim);
x0=thetac(4:end);

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