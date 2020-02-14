function f=subLhoods(NNbar,params,xdata,ydataNX)
%Input data - hospitalisations (i.e. without multipliers)

%Fit to subset of data:
threshl=17;%Lower (included in fit)
threshu=52;%Upper
ydataNX(xdata<threshl,:)=[];
xdata(xdata<threshl)=[];

ydataNX(xdata>threshu,:)=[];
xdata(xdata>threshu)=[];

totalInc=0;
lx=length(xdata);
%NNbar=[1464566;3790730;10236474;3984200;2412129];
nbar=length(NNbar);
%ydataNX=ydataNX(xdata,:);

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];

ysim=subPandemicSimulation(NNbar,params,xdata,0,0,0,243);
%{
%Eliminate age group(s):
ageOut=5;
nbar=nbar-length(ageOut);
ydataNX(:,ageOut)=[];
ysim(:,ageOut)=[];
mu(ageOut)=[];
sig(ageOut)=[];
%}
%ysim=ysim(xdata,:);%simCut in cdcPandemicSimulationW5
if totalInc==1
    muvec=sum(mu);
    sigvec=sqrt(sum(sig.^2));
    hosp=sum(ydataNX,2);
    y=sum(ysim,2);
else
    muvec=reshape(repmat(mu,lx,1),nbar*lx,1);
    sigvec=reshape(repmat(sig,lx,1),nbar*lx,1);
    hosp=reshape(ydataNX,nbar*lx,1);
    y=reshape(ysim,nbar*lx,1);  
end
%}
y=y./hosp;
y(hosp==0)=0;
%
thresh=10;
x=zeros(lx*nbar,1);
%x(y<thresh)=1;
x(hosp<thresh)=1;
y(x==1)=NaN;
if totalInc==0
    muvec(x==1)=NaN;
    sigvec(x==1)=NaN;
end
%}
%
%Normal likelihood:
L=-log(normpdf(y,muvec,sigvec));
f=nansum(L);
%}
%{
figure
L2=-log(normpdf(zeros(length(y),1),muvec,sigvec));
if totalInc==0
    L=reshape(L,lx,nbar);
    L2=reshape(L2,lx,nbar);
end
X=[L;L2];
X(isinf(X)==1)=0;
maxy=max(max(X));
fs=12; lw=2;
hold on
plot(xdata,L,'x','markersize',4,'linewidth',lw);
plot(xdata,L2,':','linewidth',lw);
xlabel('Time (weeks)','FontSize',fs);
ylabel('Likelihood','FontSize',fs);
set(gca,'FontSize',fs);
axis([1,xdata(end),0,maxy]);
if totalInc==1
    legend('L(sim)','L(0)','location','SE')
end
%legend('0-4','5-17','18-49','50-64','65+','location',NE)
grid on
grid minor
box on
hold off
%}
end