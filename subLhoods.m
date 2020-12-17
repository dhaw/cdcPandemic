function f=subLhoods(NNbar,params,xdata,ydataNX,vaxparams)
%Input data - hospitalisations (i.e. without multipliers)
totalInc=0;
%thresh=floor(max(max(ydataNX))/10);
z1combine=0;%Requires totalInc=0 & thresh=0
ages2and3=0;
mumult=1;
sigmult=1;
%Fit to subset of data:
%
threshl=17;%23;%17;%17;%Lower (included in fit)
threshu=34;%52;%34;%52;%Upper
tswitch=35-threshl+1;

ydataNX(xdata<threshl,:)=[];
xdata(xdata<threshl)=[];
ydataNX(xdata>threshu,:)=[];
xdata(xdata>threshu)=[];

maxinc=max(max(ydataNX));
thresh=1*floor(max(max(ydataNX))/10);

%{
ydataNX(11:15,:)=[];%thresh1=21 %10:16
xdata(11:15,:)=[];
%}
lx=length(xdata);
%NNbar=[1464566;3790730;10236474;3984200;2412129];
nbar=length(NNbar);
%ydataNX=ydataNX(xdata,:);

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];%*2.74;%*2.74 for under-reporting
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
%
mu=mu/mumult;
sig=sig*sigmult;
%}
%try
    %ysim=subPandemicSimulation(NNbar,params,xdata,0,0,0,243);
    ysim=subPandemicSimulationVax(NNbar,params,xdata,0,0,0,243,vaxparams);
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
    if z1combine==1 && totalInc==0
        yw1=ysim(1:tswitch,:);
        yw1=sum(yw1,1);
        muw1=mu;%*tswitch;
        sigw1=sig;%*tswitch;
        hospw1=ydataNX(1:tswitch,:);
        hospw1=sum(hospw1,1);
        %
        yw2=[yw1;ysim(tswitch+1:end,:)];
        lw2=size(yw2,1)-1;
        muvec=reshape([muw1;repmat(mu,lw2,1)],nbar*(lw2+1),1);
        sigvec=reshape([sigw1;repmat(sig,lw2,1)],nbar*(lw2+1),1);
        hospSum=repmat(sum([hospw1;ydataNX(tswitch+1:end,:)],2),nbar,1);
        hosp=reshape([hospw1;ydataNX(tswitch+1:end,:)],nbar*(lw2+1),1);
        y=reshape(yw2,nbar*(lw2+1),1); 
    elseif totalInc==1
        muvec=sum(mu);%Weight by populations??
        sigvec=sqrt(sum(sig.^2));
        hosp=sum(ydataNX,2);
        y=sum(ysim,2);
    else
        muvec=reshape(repmat(mu,lx,1),nbar*lx,1);
        sigvec=repmat(sig,lx,1);

        %sigvec(9:15,:)=sigvec(9:15,:)*10;
        %sigvec(sigvec<.5*maxinc)=2*sigvec(sigvec<.5*maxinc);

        sigvec=reshape(sigvec,nbar*lx,1);
        hospSum=repmat(sum(ydataNX,2),nbar,1);
        hosp=reshape(ydataNX,nbar*lx,1);
        y=reshape(ysim,nbar*lx,1);
        
        xrep=repmat(xdata,nbar,1);
        %sigvec(hosp<30)=10*sigvec(hosp<30);
    end%
    %}
    y=y./hosp;
    y(hosp==0)=0;
    %
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
    toLog=normpdf(y,muvec,sigvec);
    if ages2and3==1
        toLog=reshape(toLog,lx,nbar);
        toLog(xdata<35,[1,4,5])=NaN;%Assumes na=5
        toLog=reshape(toLog,nbar*lx,1);
    end
    L=-log(toLog);
    f=nansum(L);
    %}
%catch
    %f=Inf;
%end
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