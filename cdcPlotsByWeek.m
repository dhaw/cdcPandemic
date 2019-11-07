function f=cdcPlotsByWeek(x)
%Needs weekly data!

%x=sortrows(x,[2,3]);

inds09=find(x(:,2)==2009);
inds10=find(x(:,2)==2010);
inds11=find(x(:,2)==2011);
x=x([inds09;inds10;inds11],:);
x=sortrows(x,[2,3,4]);
years=x(1,2):x(end,2);
%numyears=length(years);

%
twelve=kron((1:12)',ones(4,1));%No 5th weeks in data
four=1:4;
numyears=length(years);
y=[kron(years',ones(48,1)),repmat(twelve,numyears,1),repmat(four',12*numyears,1)];
y(:,4:18)=nan;
%Clunky way:

for i=1:numyears
    xi=x(x(:,2)==years(i),:);
    for j=1:12
        %month=find(xi(:,3)==j);
        xij=xi(xi(:,3)==j,:);
        if isempty(xij)==0
            xij=sortrows(xij,4);
            yind=(i-1)*48+(j-1)*4+1;
            %xij=nanmean(xi(month,:),1);
            y(yind:yind+size(xij,1)-1,4:end)=xij(:,5:end);
        end
    end
end
%}
%y=x(:,3:end);
incScaled=y(:,7)./y(:,8).*y(:,17);
tvec=1:size(incScaled,1);

tvec1=y(:,2); tvec2=y(:,3);
%tvec(isnan(incScaled))=[];
%tvec1(isnan(incScaled))=[];
%tvec2(isnan(incScaled))=[];
%incScaled(isnan(incScaled))=[];
jans=find(tvec1==1);
jans1=find(tvec2==1);
jans=intersect(jans,jans1);
jans=tvec(jans);
ljans=length(jans);

fs=12; lw=2;
ticks=jans;
labs={'01/09','01/10','01/11'};

figure
hold on
mx=max(incScaled);
%
for j=1:ljans
    plot([jans(j),jans(j)],[0,mx],'k--','linewidth',lw)
end
%}
plot(tvec,incScaled,'-','linewidth',lw)
%Manual (for plotting):
plot(tvec(19:48),incScaled(19:48),'-','linewidth',lw+1)
hold off
xlabel('Time (weeks)')
ylabel('Incidence')
%xvec=1:6:tvec(end);
set(gca,'xtick',ticks,'xticklabel',labs)
set(gca,'fontsize',fs)
ax=gca;
ax.XAxis.MinorTickValues=tvec(1:4:end);
axis([0,112,0,mx])%tvec(end)
grid on
grid minor
box on

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