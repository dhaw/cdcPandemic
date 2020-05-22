function f=plotz2(cellZ2,real)
numplots=8;
%cell in alphabetical order of 2-letter state code
if numplots==8
    names={'CA','CO','CT','GA','MD','NY','OR','TN'};%Not MN/NM/US
elseif numplots==6
    names={'CO','CT','GA','MD','NY','OR','TN'};%Not CA/MN/NM/US
end
%,'NM' - no xnm
hospUnder=2.74;
hospUnderMin=1.7;
hospUnderMax=4.5;

legString={'0-4','5-17','18-49','50-64','65+'};
cmap=lines(7);
transp=.1;
fs=10; lw=2;
figure('DefaultAxesPosition', [0.1, 0.1, 0.875, 0.875])
%sgtitle('Hospitalisations/1000 (prior: USA)','fontsize',fs)
for i=1:numplots
    h=subplot(numplots/2,2,i);
    ci=cellZ2{i};
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
    if i>numplots-2
        xlabel('Age group')
    end
    if rem(i,2)==1%i==1 || i==numplots/2+1
        ylabel('Hospitalisations')
    else
        yticks(0:.5:1)
        yticklabels({'',''})
    end
    title(names{i})
    set(gca,'fontsize',fs)
    xticks(1:nbar)
    xticklabels(legString)
    axis([0,nbar+1,0,1.2])
    %set(h,'ActivePositionProperty','OuterPosition')%,[0.00    0.75    0.33    0.25])
    %legend([hleg1,hleg2,hleg3,hleg4,hleg5],legString,'location','NE')[
    grid on
    grid minor
    box on
end