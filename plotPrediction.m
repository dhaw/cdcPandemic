function f=plotPrediction(pd)
NNbar=[19169690;62121035;184015269;39570590];
tswitch=243;

pvals=[1.4];
[~,outvec]=pandemic1DallV(pvals,4:7,0,0,0);%Output tvec and incidence
%Absolute incidence
tvec=outvec(:,1); ivec=outvec(:,2);
ivec(tvec>tswitch)=[];
tvec(tvec>tswitch)=[];
plot1=[tvec,ivec];

burn=200;
pd=pd(burn+1:end,:);
%paramsMLE=[];
vals=prctile(pd,[5,25,50,75,95]);
lv=length(vals);
c=cell(lv,1);
maxi=max(ivec);
for i=1:lv
    [~,outvec]=pandemic1DallV(vals(i),4:7,0,0,0);%Output tvec and incidence
    %Absolute incidence
    tvec=outvec(:,1); ivec=outvec(:,2);
    ivec(tvec<tswitch)=[];
    tvec(tvec<tswitch)=[];
    c{i}=[tvec,ivec];
    maxi=max(maxi,max(ivec));
end
fs=12; lw=2;
cmap=repmat([.8,.5,.2,.5,.8]',1,3);
figure;
hold on
plot(plot1(:,1),plot1(:,2),'k-','linewidth',lw)
for i=1:lv
    ploti=c{i};
    plot(ploti(:,1),ploti(:,2),'-','linewidth',lw,'color',cmap(i,:))
end
xlabel('Time (days)')
ylabel('Incidence')
set(gca,'fontsize',fs)
axis([0,400,0,maxi])
grid on
grid minor
box on