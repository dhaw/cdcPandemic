function f=cdcPlotsAge(x)
%x=sortrows(x,[2,3]);
byAge=1;

inds09=find(x(:,2)==2009);
inds10=find(x(:,2)==2010);
inds11=find(x(:,2)==2011);
x=x([inds09;inds10;inds11],:);
years=x(1,2):x(end,2);
%numyears=length(years);

twelve=kron((1:12)',ones(4,1));%No 5th weeks in data
four=1:4;
numyears=length(years);
y=[kron(years',ones(48,1)),repmat(twelve,numyears,1),repmat(four',12*numyears,1)];
y(:,4:18)=nan;
%y(:,3:17)=nan;
%Clunky way:
for i=1:numyears
    xi=x(x(:,2)==years(i),:);
    for j=1:12
        month=find(xi(:,3)==j);
        if isempty(month)==0%Assumes all for age groups present for each data point
            xij=xi(month,:);
            yind=(i-1)*48+(j-1)*4;
            for k=1:4
                age=xij(xij(:,4)==k,:);
                age=nanmean(age,1);
                y(yind+k,4:end)=age(5:end);
            end
        end
    end
end

if byAge==1
    incScaled=y(:,7)./y(:,8).*y(:,17);%.*y(:,15);%.*y(:,6);%.*y(:,end-1);
    tvec1=y(1:4:end,1); tvec2=y(1:4:end,2);
    tvec=1:size(incScaled,1)/4;

    %tvec(isnan(incScaled))=[];
    %tvec1(isnan(incScaled))=[];
    %tvec2(isnan(incScaled))=[];
    %incScaled(isnan(incScaled))=[];

    inc1=incScaled(1:4:end);
    inc2=incScaled(2:4:end);
    inc3=incScaled(3:4:end);
    inc4=incScaled(4:4:end);

    f=[inc1,inc2,inc3,inc4];

    jans=tvec(tvec2==1);%-.5;
    ljans=length(jans);
    fs=12; lw=2;
    ticks=[1,13,25];
    labs={'01/09','01/10','01/11'};
    cmap=lines(4);
    %
    figure
    hold on
    mx=max(incScaled);
    for j=1:ljans
        plot([jans(j),jans(j)],[0,mx],'k--','linewidth',lw)
    end
    h1=plot(tvec,inc1,'-','linewidth',lw,'color',cmap(1,:));
    h2=plot(tvec,inc2,'-','linewidth',lw,'color',cmap(2,:));
    h3=plot(tvec,inc3,'-','linewidth',lw,'color',cmap(3,:));
    h4=plot(tvec,inc4,'-','linewidth',lw,'color',cmap(4,:));
    %Manual (for plotting):
    %plot(tvec(5:12),incScaled(5:12),'-','linewidth',lw+1)
    hold off
    xlabel('Time (months)')
    ylabel('Incidence (relative)')
    %ylabel('Incidence (absolute)')
    xvec=1:6:tvec(end);
    set(gca,'xtick',ticks,'xticklabel',labs)
    set(gca,'fontsize',fs)
    ax=gca;
    ax.XAxis.MinorTickValues=tvec;
    axis([0,28,0,mx])%tvec(end)
    legend([h1,h2,h3,h4],'0-4','5-19','20-64','65+')
    grid on
    grid minor
    box on
else
    incScaled=y(:,7).*y(:,17);%.*y(:,15);%.*y(:,6);%.*y(:,end-1);
    pop=y(:,8);
    pop1=pop(1:4:end);
    pop2=pop(2:4:end);
    pop3=pop(3:4:end);
    pop4=pop(4:4:end);
    pop=pop1+pop2+pop3+pop4;
    
    tvec1=y(1:4:end,1); tvec2=y(1:4:end,2);
    tvec=1:size(incScaled,1)/4;

    %tvec(isnan(incScaled))=[];
    %tvec1(isnan(incScaled))=[];
    %tvec2(isnan(incScaled))=[];
    %incScaled(isnan(incScaled))=[];

    inc1=incScaled(1:4:end);
    inc2=incScaled(2:4:end);
    inc3=incScaled(3:4:end);
    inc4=incScaled(4:4:end);
    inc=(inc1+inc2+inc3+inc4)./pop;

    f=inc;

    jans=tvec(tvec2==1);%-.5;
    ljans=length(jans);
    fs=12; lw=2;
    ticks=[1,13,25];
    labs={'01/09','01/10','01/11'};
    %
    figure
    hold on
    mx=max(inc);
    for j=1:ljans
        plot([jans(j),jans(j)],[0,mx],'k--','linewidth',lw)
    end
    plot(tvec,inc,'-','linewidth',lw);
    %Manual (for plotting):
    %plot(tvec(5:12),incScaled(5:12),'-','linewidth',lw+1)
    hold off
    xlabel('Time (months)')
    ylabel('Incidence (relative)')
    %ylabel('Incidence (absolute)')
    set(gca,'xtick',ticks,'xticklabel',labs)
    set(gca,'fontsize',fs)
    ax=gca;
    ax.XAxis.MinorTickValues=tvec;
    axis([0,14,0,mx])%tvec(end)
    grid on
    grid minor
    box on
end
%{
twelve=1:12;
years=x(1,2):x(end,2);
numyears=length(years);
months=[kron(years',ones(12,1)),repmat(twelve',numyears,1)];
weeks=[kron(months,ones(4,1)),repmat((1:4)',length(months),1)];

dummy=100*weeks(:,1)+weeks(:,2)+weeks(:,3)/4-.25;
y=weeks; y(:,4:size(x,2)-1)=nan;
%Clunky:
xyear=x(:,2); xmonth=x(:,3); xweek=x(:,4);
xdummy=100*xyear+xmonth+xweek/4-.25;
for i=1:length(xdummy)
    %ind=find(dummy==xdummy(i));
    y(dummy==xdummy(i),4:end)=x(i,5:end);
end

tvec2=num2str(tvec2); tvec2(isspace(tvec2))='0';
tlabs=strcat(num2str(tvec1),'/',tvec2);
tlabs=tlabs(:,3:end);
%tlabs=tlabs(find(~isspace(tlabs)));
%}