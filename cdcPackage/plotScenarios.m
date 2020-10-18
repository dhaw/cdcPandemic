function f=plotScenarios(F,NNbar)
numScen=4;
NN=sum(NNbar);
%
fi=F{numScen};
f1=median(sum(fi(:,1:5,:),2),3)/NN*1e5;
f2=median(sum(fi(:,6:10,:),2),3)/NN*1e5;
f3=median(sum(fi(:,11:15,:),2),3)/NN*1e5;
f4=median(sum(fi(:,16:20,:),2),3)/NN*1e5;
lf=size(f1,1);
%}
%{
f1=median(sum(F{numScen}(:,1:5,:),2),3)/NN*1e5;
f2=median(sum(F{numScen+1},2),3)/NN*1e5;
f3=median(sum(F{numScen+2},2),3)/NN*1e5;
f4=median(sum(F{numScen+3},2),3)/NN*1e5;
lf=size(f1,1);
%}
fs=10; lw=2;
cmap=lines(7);
figure
hold on
h1=plot(1:lf,f1,'linewidth',lw,'color',cmap(1,:));
h2=plot(1:lf,f2,'-','linewidth',lw,'col',cmap(2,:));
h3=plot(1:lf,f3,'-','linewidth',lw,'col',cmap(3,:));
h4=plot(1:lf,f4,'-','linewidth',lw,'col',cmap(4,:));
legend([h1,h2,h3,h4],'SI','HO','DE','AV');
%%
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel('Incidence')
axis tight%([1,52,0,6.6e6])
grid on
grid minor
box on
f=0;
%%
%{
function f=plotScenarios(F)
numScen=14;
%
ivec=1:numScen;
z=zeros(1,numScen);
f1=F{ivec(1)};
f1=sum(f1(:,1:5,1),2);
lf=length(f1);
X=zeros(lf,numScen);
X(:,1)=f1;
z(1)=sum(f1);
h=zeros(1,6);
for i=2:numScen
    fi=F{ivec(i)};
    X(:,i)=sum(fi(:,1:5,1),2);
    z(i)=sum(sum(fi(:,:,1)));
end
f=z;
%
fs=10; lw=2;
cmap=repmat(lines(6),4,1);
figure
hold on
%plot(1:lf,X,'linewidth',lw);
%
cmap=lines(2);
for i=1:numScen
    %h(i)=plot(1:lf,X(:,i),'linewidth',lw,'color',cmap(i,:));
    if i<=numScen/2
        col=cmap(1,:);
    else
        col=cmap(2,:);
    end
    plot(1:lf,X(:,i),'linewidth',lw,'color',col);%,cmap(i,:));
end
h1=plot(-1,-1,'-','linewidth',lw,'col',cmap(1,:));
h2=plot(-1,-1,'-','linewidth',lw,'col',cmap(2,:));
legend([h1,h2],'EPI1','EPI2','location','NE');
maxx=max(max(X));
%}
%%
%{
fs=10; lw=2;
cmap=lines(6);
figure
%plot(1:lf,X,'linewidth',lw);
%
lf=size(F{1},1);
tvec=1:lf;
tvec2=[tvec,fliplr(tvec)];
maxx=1;
ivec=[1,8]+4;
h=zeros(1,length(ivec));
hold on
for i=1:length(ivec)%numScen
    X=F{ivec(i)};
    X=prctile(X,[25,50,75],2);
    y1=X(:,1);
    y2=X(:,2);
    y3=X(:,3);
    maxy=max(max(X));
    %if ismember(ivec(i),[4,6,8,11,16,18,21,23])
        %plot([0,53],.01*21888099*[1,1],'k--','linewidth',1)%Sum(NNbarUS)
    %end
    %if ismember(ivec(i),[4,6,8,11,16,18,21,23]+1)
        %plot(34+8.7143*[1,1],[0,maxy],'k--','linewidth',1)
    %end
    %if ismember(ivec(i),[2,4,5,9,10,14,21,22])
        %plot((243+107*[1,1])/7,[0,6.6e6],'k--','linewidth',1)
    %end
    inBetween=[y1;flipud(y3)];
    fill(tvec2,inBetween,cmap(i,:),'facealpha',.2);
    plot(tvec,y1,'-','linewidth',1,'color',cmap(i,:));
    plot(tvec,y3,'-','linewidth',1,'color',cmap(i,:));
    h(i)=plot(tvec,y2,'-','linewidth',2,'color',cmap(i,:));
    maxx=max(maxx,maxy);
end
%}
%%
%{
set(gca,'fontsize',fs)
xlabel('Time (weeks)')
ylabel('Incidence')
axis([1,52,0,6.6e6])%maxx])
%legend(h,num2str((ivec)'))
grid on
grid minor
box on
%}
f=0;