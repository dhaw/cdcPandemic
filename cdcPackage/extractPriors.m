function [f,g]=extractPriors(xsto,x0)
burn=1000;
%
%Gamma:
lx=size(xsto,2);
priors=zeros(2,lx);
for i=1:lx
    pd=fitdist(xsto(burn+1:end,i),'gamma');
    priors(1,i)=pd.a;
    priors(2,i)=pd.b;
end
f=priors;
g=[];
%}
%{
%Histogram:
lx=size(xsto,1);
plim=[1,0;
1,0;
1,0;
240,0;
2,1;%max,min
2.4,.6]';
Cc=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
Cvec=reshape(x0(3:27),25,1);
%Cvec=thetac(3:27)';
plim=[[10*Cvec,.1*Cvec];plim'];
lp=length(plim);

pdiff=plim(:,1)-plim(:,2);
steps=0:.01:1;
ls=length(steps);
prows=repmat(plim(:,2),1,ls)+kron(steps,pdiff);

priors=zeros(lp,ls-1);
binEdges=zeros(lp,ls);
for i=1:lp
    pd=histogram(xsto(burn+1:lx,i),'binedges',prows(i,:),'normalization','pdf');
    priors(i,:)=pd.Values;
    binEdges(i,:)=pd.BinEdges;
end
f=priors;
g=binEdges;
%}