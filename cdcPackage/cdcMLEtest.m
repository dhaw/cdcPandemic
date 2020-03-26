function f=cdcMLEtest(ydataNX,thetaIn)
xdata=1:36;
lx=length(xdata);
mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
nbar=length(mu);
muvec=reshape(repmat(mu,lx,1),nbar*lx,1);
sigvec=reshape(repmat(sig,lx,1),nbar*lx,1);

nbar=size(ydataNX,2);
ydataNX=ydataNX(xdata,:);
hosp=reshape(ydataNX,nbar*lx,1);

thresh=100;
x=zeros(lx*nbar,1);
x(hosp<thresh)=1;

options=optimset('MaxFunEvals',100000);%,'TolX',10^-30,'TolFun',10^-30);
%[theta,~]=mle(hosp,'logpdf',@(hosp,params)lhoodFct1(hosp,params,xdata,lx,nbar,x,muvec,sigvec),'start',thetaIn);
theta=fminsearch(@(params)lhoodFct1(hosp,params,xdata,lx,nbar,x,muvec,sigvec),thetaIn,options);

f=theta;
end

function f=lhoodFct1(hosp,params,xdata,lx,nbar,x,muvec,sigvec)
%{
xdata=1:36;
lx=length(xdata);
mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
nbar=length(mu);
muvec=reshape(repmat(mu,lx,1),nbar*lx,1);
sigvec=reshape(repmat(sig,lx,1),nbar*lx,1);
nbar=size(ydataNX,2);
thresh=100;
x=zeros(lx*nbar,1);
x(hosp<thresh)=1;
%}
ysim=cdcPandemicSimulationW5(params,xdata,0,0,0,243);
y=reshape(ysim,nbar*lx,1);  
y(hosp==0)=0;
y=y./hosp;
y(x==1)=NaN;
%muvec(x==1)=NaN;
%sigvec(x==1)=NaN;
f=-nansum(log(normpdf(y,muvec,sigvec)));
end