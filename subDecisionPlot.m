function f=subDecisionPlot(z)

x=(1:.05:2);%xi - match cdcDecisionSeverity
lx=length(x);
points=size(z,2);

figure
fs=10; lw=2; ms=5;
transp=.1;
xx=ones(points,1);
cmap=lines(7);
maxy=max(max(max(z)));
NstarVec=200:200:16000;%12000;
ln=length(NstarVec);

%
%for n=1:ln
    Nstar=1600;%NstarVec(n);
    
    zn=z;
    for i=1:lx
        zn(i,zn(i,:,1)>Nstar,1)=NaN;
        zn(i,zn(i,:,1)<=Nstar,2)=NaN;
    end
    %figure
    %subplot(2,5,n)
    hold on
    plot(x,Nstar*ones(1,lx),'k--','linewidth',lw)
    for i=1:lx
        h1=scatter(x(i)*xx,zn(i,:,1),'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:));
        h1.MarkerFaceAlpha=transp;
        h1.MarkerEdgeAlpha=transp;
        h2=scatter(x(i)*xx,zn(i,:,2),'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:));
        h2.MarkerFaceAlpha=transp;
        h2.MarkerEdgeAlpha=transp;
    end
    xlabel('Relative severity')
    ylabel('Hospitalisations')
    title(strcat('N^*=',num2str(Nstar)))
    set(gca,'fontsize',fs)
    axis([x(1),x(end),0,maxy])
    grid on
    grid minor
    box on
%end
%}

%{
NstarVec=800:10:1400;
ln=length(NstarVec);

[a,b,c]=size(z);
zn=zeros(b,c);
zn(:,:)=z(1,:,:);%Select severity
zplot=repmat(zn(:,1),1,ln);
zexp=mean(zn(:,1));
    for n=1:ln
        if NstarVec(n)>zexp
            zplot(:,n)=zn(:,2);
        end
    end
    hold on
    for n=1:ln
        h1=scatter(NstarVec(n)*ones(b,1),zplot(:,n),'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:));
        h1.MarkerFaceAlpha=transp;
        h1.MarkerEdgeAlpha=transp;
        %h2=scatter(x(i)*xx,zn(i,:,2),'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:));
        %h2.MarkerFaceAlpha=transp;
        %h2.MarkerEdgeAlpha=transp;
    end
    xlabel('N^*')
    ylabel('Hospitalisations')
    %title(strcat('N^*=',num2str(Nstar)))
    set(gca,'fontsize',fs)
    axis([NstarVec(1),NstarVec(end),0,maxy])
    grid on
    grid minor
    box on
%}
    
%{
figure
fs=12; lw=2; ms=3;
transp=.1;
cmap=lines(7);
maxy=max(max([Cfsd;Cgsd]));
hold on
plot([35,35],[0,maxy],'k--','linewidth',lw)
for i=1:size(Cfsd,2)
    h1=scatter(param',Cfsd(:,i),'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:));
    h1.MarkerFaceAlpha=transp;
    h1.MarkerEdgeAlpha=transp;
end
for i=1:size(Cgsd,2)
    h2=scatter(param',Cgsd(:,i),'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:));
    h2.MarkerFaceAlpha=transp;
    h2.MarkerEdgeAlpha=transp;
end
xlabel('\epsilon')%X
%ylabel('cc_{pand}')
ylabel('\sigma(cc_{pand})')
set(gca,'fontsize',fs)
axis([param(1),param(end),0,maxy])%X
grid on
grid minor
box on
%}