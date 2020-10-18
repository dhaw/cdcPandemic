function [pr,NN,n,nbar,na,NNbar,NNrep,Dout,beta]=hePrepCovid19(D)%,hvec,muvec)
%D - column vector or populations.
%Possible generalsiation to within-sector heterogeneity - one column per
%category. 
%hvec - hospitalisations by age group. 
%muvec - deaths by age group. 
%Currently 1 age group, so hvec and muvec are scalars. 

urbrur=0;%Turn in to @home vs @work
%%
%Population density:
[n,na]=size(D);
nbar=n*na;
NN=sum(D,2);
NNbar=reshape(D,nbar,1);
NNrep=repmat(NN,na,1);
%{
if urbrur==1
    hvec=kron(hvec,ones(n,1));
    muvec=kron(muvec,ones(n,1));
    %Kkron=[1,.05;.05,.75];%Sonoma
    Kkron=[1,.05;.05,.85];%Mendocino
else
    Kkron=1;
end
%}
%{
C=[.3827    0.9115    0.0419;
    1.2062    7.3538    0.1235;
    0.1459    0.4810    0.1860];
%}
C=eye(na);

%K=heMakeDs(NN,eye(10));
K=heMakeDs(NN,ones(length(NN)-1,1));%,1);

%K=rand(n);
%K=normr(K);
D=kron(C,K);
%%
%toHosp=3;%Symp to hosp%****
propHospGivenInf=.04;
Thosp=5;%In hosp****
Text=4.6;
Tonset=1;
pr=struct;
pr.sigma=1/Text;
pr.omega=1/Tonset;
pr.g1=1/2.1;
pr.g2=1/2.1;
pr.g3=1/Thosp;%****
pr.gX=1/4;
pr.p1=1-.34;
pr.p2=.505;%.04; hvec;
pr.p3=1;
pr.p4=1;
pr.q1=0;
pr.q2=0;
%
pr.g4=1/(1/pr.g2-1/pr.q1);%~q's and other gammaspr.p1=1-.34;
pr.g4X=1/(1/pr.gX+1/pr.q2);
%
ph=propHospGivenInf/pr.p1/pr.p2;%proportion of severe that are hospitalised%****
pr.h=ph/(1-ph);%1/toHosp;
pr.mu=.0025/.04;%muvec;
pr.odds=0;
pr.qnew=0;
pr.red=2/3;
pr.R0=2.5;
%%
%Debug stuff:
%hvec=zeros(nbar,1);
%gammaEff=(1-pr.p1)*pr.g1*ones(na*n,1)+pr.p1*((1-pr.h)*pr.g2+hvec*pr.p2);
%Vinv=diag(1./gammaEff);
%
%%
%isdual=1 - equivalent here
%Ceff=kron(C,Ckron);%Urb/rural mixing here
Dout=D;
Deff=Dout.*repmat(NNbar,1,n*na)./repmat(NNbar',n*na,1);
%
ntot=n*na;
onesn=ones(ntot,1);
F=zeros(4*ntot,4*ntot);
%F(1:ntot,1:ntot)=Deff;
F(1:ntot,ntot+1:end)=[pr.red*Deff,repmat(Deff,1,2)];
vvec=kron([pr.sigma;pr.g1;pr.g2;pr.gX+pr.h],ones(ntot,1));
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*(1-pr.p2)*onesn);
V(3*ntot+1:4*ntot,1:ntot)=diag(-pr.sigma*pr.p1*pr.p2*onesn);
GD=F/V;
%{
%HE:
ntot=n*na;
F=zeros(7*ntot,7*ntot);
F(1:ntot,1:ntot)=Deff;%ntot+1:end)=[repmat(2/3*Deff,1,2),repmat(Deff,1,4)];
onesn=ones(ntot,1);
vvec=[pr.sigma*onesn;pr.g1*onesn;pr.omega*onesn;pr.g2*onesn;pr.g2*onesn;pr.h*onesn;pr.h*onesn];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*onesn);
V(3*ntot+1:4*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p3).*(1-pr.p2));
V(4*ntot+1:5*ntot,2*ntot+1:3*ntot)=diag(-pr.p3.*(1-pr.p2));
V(5*ntot+1:6*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p4).*pr.p2);
V(6*ntot+1:7*ntot,2*ntot+1:3*ntot)=diag(-pr.p4.*pr.p2);
GD=F/V;
%}
d=eigs(GD,1); R0a=max(d); beta=pr.R0/R0a;
end