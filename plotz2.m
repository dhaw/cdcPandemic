function f=plotz2(cellZ2,real)
numplots=10;
%cell in alphabetical order of 2-letter state code
if numplots==8
    names={'CA','CO','CT','GA','MD','NY','OR','TN'};%Not MN/NM/US
elseif numplots==6
    names={'CO','CT','GA','MD','NY','OR','TN'};%Not CA/MN/NM/US
elseif numplots==10
    names={'CA','CO','CT','GA','MD','MN','NM','NY','OR','TN'};%Not CA/MN/NM/US
end
%,'NM' - no xnm
hospUnder=1;%2.74;
hospUnderMin=1.7/2.74;
hospUnderMax=4.5/2.74;

%legString={'0-4','5-17','18-49','50-64','65+'};
legString={'0-4','-17','-49','-64','65+'};
cmap=lines(7);
transp=.1;
fs=10; lw=2;
figure('DefaultAxesPosition', [0.1, 0.1, 0.875, 0.875])
%sgtitle('Hospitalisations/1000 (prior: USA)','fontsize',fs)
for i=1:numplots
    h=subplot(2,numplots/2,i);
    ci=cellZ2{i};
    %
    %Comment out for Z1
    ci=ci(10:5:end,:)';
    mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
    sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
    repmu=repmat(mu',1,size(ci,2));
    repsig=repmat(sig',1,size(ci,2));
    ci=ci./normrnd(repmu,repsig);
    %}
    reali=real(:,i);
    maxy=max(max(ci));
    maxy=max(maxy,max(reali));
    [nbar,lp]=size(ci);
    onesi=ones(lp,1);
    hold on
    violinplot(ci');
    %plot(1:nbar,reali,'kx','markersize',10,'linewidth',lw)
    %
    for j=1:nbar
        %h1=scatter(j*onesi,ci(j,:),'markerfacecolor',cmap(j,:),'markeredgecolor',cmap(j,:));
        %h1.MarkerFaceAlpha=transp;
        %h1.MarkerEdgeAlpha=transp;
        %plot(j,reali(j),'kx','markersize',10,'linewidth',lw)%'markerfacecolor',cmap(j,:),'markeredgecolor',cmap(j,:));
        plot(j+.1,hospUnder*reali(j),'kx','markersize',10,'linewidth',lw)
        plot([j,j]+.1,reali(j)*[hospUnderMin,hospUnderMax],'k-','linewidth',lw)
    end
    %}
    %{
    hleg1=scatter(-1,-1,'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:));
    hleg2=scatter(-2,-2,'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:));
    hleg3=scatter(-3,-3,'markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:));
    hleg4=scatter(-4,-4,'markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:));
    hleg5=scatter(-5,-5,'markerfacecolor',cmap(5,:),'markeredgecolor',cmap(5,:));
    %}
    hold off
    if i>numplots/2%/2%-2
        xlabel('Age group')
    end
    if i==1 || i==numplots/2+1%i==1 || i==numplots/2+1%rem(i,2)==1%i==1 || i==numplots/2+1
        ylabel('Hosp./1000')
    else
        yticks(0:.5:1)
        yticklabels({'',''})
    end
    title(names{i})
    set(gca,'fontsize',fs)
    xticks(1:nbar)
    xticklabels(legString)
    axis ([0,nbar+1,0,1])
    %set(h,'ActivePositionProperty','OuterPosition')%,[0.00    0.75    0.33    0.25])
    %legend([hleg1,hleg2,hleg3,hleg4,hleg5],legString,'location','NE')[
    grid on
    grid minor
    box on
end