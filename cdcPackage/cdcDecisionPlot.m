function f=cdcDecisionPlot(z)
x=(1:.05:2);%xi - match cdcDecisionSeverity
lx=length(x);
points=size(z,2);

figure
fs=10; lw=1; ms=5;
xx=ones(points,1);
cmap=lines(7);
maxy=max(max(max(z)));
NstarVec=6000:1000:15000;%12000;
ln=length(NstarVec);
for n=1:ln
    Nstar=NstarVec(n);
    zn=z;
    for i=1:lx
        zn(i,zn(i,:,2)>Nstar,1)=NaN;
        zn(i,zn(i,:,2)<=Nstar,2)=NaN;
    end
    %{
    fs=12; lw=2; ms=5;
    xx=ones(points,1);
    cmap=lines(7);
    maxy=max(max(max(z)));
    %}
    %figure
    subplot(2,5,n)
    hold on
    for i=1:lx
        plot(x(i)*xx,zn(i,:,1),'o','color',cmap(1,:),'markersize',ms,'linewidth',lw)
        plot(x(i)*xx,zn(i,:,2),'o','color',cmap(2,:),'markersize',ms,'linewidth',lw)
    end
    xlabel('Relative severity')
    ylabel('Hospotalisations')
    title(strcat('N^*=',num2str(Nstar)))
    set(gca,'fontsize',fs)
    axis([x(1),x(end),0,maxy])
    grid on
    grid minor
    box on
end

%{
for n=1:length(Nvec)
    
end
%}