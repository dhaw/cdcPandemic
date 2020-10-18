function f=plotDataSub(xx,yyx,NN)
numplots=length(xx);
names={'CA','CO','CT','GA','MD','MN','NM','NY','OR','TN'};
fs=10; lw=2;
figure('DefaultAxesPosition', [0.1, 0.1, 0.875, 0.875])
for i=1:numplots
    NNi=sum(NN{i});
    xi=xx{i};
    yi=yyx{i}/NNi*10000;
    h=subplot(numplots/2,2,i);
    %
    %colormap(cmap5)
    plot(xi,yi,'linewidth',lw);
    %{
    hold on
    for j=1:5
        plot(xi,yi(:,j),'linewidth',lw,'color',cmap5(j,:));
    end
    hold off
    %}
    if i>numplots-2%/2%-2
        xlabel('Week')
    end
    if rem(i,2)==1%i==1 || i==numplots/2+1%rem(i,2)==1%i==1 || i==numplots/2+1
        ylabel('Hospitalisations')
    else
        yticks(0:.5:1)
        yticklabels({'',''})
    end
    title(names{i})
    set(gca,'fontsize',fs)
    %xticks(17:52)
    %xticklabels(legString)
    axis ([17,52,0,.6])
    %set(h,'ActivePositionProperty','OuterPosition')
    grid on
    grid minor
    box on
end