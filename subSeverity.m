function f=subSeverity(NNbar,xsto,xdata,ydata,vaxparams)
burn=10000;
int=500;
na=length(NNbar);
pstar=0:.05:1;
Hvec=0:5e2:1e5;
numIter=100;

lp=length(pstar);
lh=length(Hvec);


%Change for different severity:
meanHosp=[393.0256  999.3054  702.6867  406.0680  177.1958];
sdHosp=[129.6984 329.7708 231.8866 134.0024 58.4746];

M=[meanHosp;meanHosp/2;meanHosp/2];
S=[sdHosp;sdHosp/2;sdHosp/2];
tswitch=243+[0,0,28];
lm=size(M,1);
X=zeros(lm,lh);

xsto=xsto(burn+1:int:end,:);
lx=size(xsto,1);
for m=1:lm
    meanHosp=M(m,:);
    sdHosp=S(m,:);
    hospVals=zeros(lx,numIter);
    for i=1:lx
        [~,~,z2]=subPandemicSimulationVax(NNbar,xsto(i,:),xdata,0,0,ydata,tswitch(m),vaxparams);
        for j=1:numIter
            mult=normrnd(meanHosp,sdHosp);
            hospVals(i,j)=sum(z2./mult);
        end
    end
    hospVals=reshape(hospVals,lx*numIter,1);


    for h=1:lh
        X(m,h)=length(hospVals(hospVals>Hvec(h)))/lx/numIter;
    end
    
end
f=X;

fs=10; lw=2; ms=7;
cmap=lines(7);
figure;
hold on
h1=plot(Hvec,X(1,:),'linewidth',lw,'color',cmap(1,:));
h2=plot(Hvec,X(2,:),'linewidth',lw,'color',cmap(2,:));
h3=plot(Hvec,X(3,:),'--','linewidth',lw,'color',cmap(2,:));
plot([Hvec(1),Hvec(end)],[.1,.1],'k--','linewidth',1.5)
plot([7e4],[.1],'ko','markersize',ms,'markerfacecolor','k')
set(gca,'fontsize',fs)
xlabel('Hospitalisations (h)')
ylabel('P(H>h)')
axis([Hvec(1),Hvec(end),0,1])
legend([h1,h2,h3],'2009','2x severity','Delay','location','NE')
grid on
grid minor
box on