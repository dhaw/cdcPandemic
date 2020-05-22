function f=vaxTestPlot(NNbar,xdata,xsto,vaxparams)
mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
na=length(mu);
tswitchvec=(0:20);%Weeks
lt=length(tswitchvec);
%
burn=10000;
int=100;
NNtot=sum(NNbar);
prior=xsto(burn+1:int:end,:);
lp=size(prior,1);

Z=zeros(na,lt,lp);
for i=1:lt
    tswitch=7*tswitchvec(i)+243;%Days
    Zi=zeros(na,lp);
    parfor j=1:lp
        [~,~,z2]=subPandemicSimulationVax(NNbar,prior(j,:),xdata,0,0,0,tswitch,vaxparams);
        %Zi(:,j)=z2'./mu'/NNtot*1000;
        ri=normrnd(mu,sig);
        Zi(:,j)=z2'./ri'/NNtot*1000;
    end
    Z(:,i,:)=Zi;
end
%}
f=Z;
%
pandsq=prctile(Z,[25,50,75],3);
y1=pandsq(:,:,1);
y2=pandsq(:,:,2);
y3=pandsq(:,:,3);
figure
fs=10; lw=2;
cmap=lines(na);
tvec=tswitchvec;
tvec2=[tvec,fliplr(tvec)];
inBetween=[y1,fliplr(y3)];
%plot(tswitchvec,Z,'-','linewidth',lw)%,'color',col1);
h=zeros(1,na);
hold on
for i=1:na
    fill(tvec2,inBetween(i,:),cmap(i,:),'facealpha',.2);
    plot(tvec,y1(i,:),'-','linewidth',1,'color',cmap(i,:));
    plot(tvec,y3(i,:),'-','linewidth',1,'color',cmap(i,:));
    h(i)=plot(tvec,y2(i,:),'-','linewidth',2,'color',cmap(i,:));
end
xlabel('Delay (weeks)','FontSize',fs);
ylabel('h_2/1000 population');
set(gca,'FontSize',fs);
maxY=max(max(y3));
axis ([0,lt-1,0,maxY])
legend(h,{'0-4','5-17','18-49','50-64','65+'},'location','NE')
grid on
grid minor
box on
hold off
%}