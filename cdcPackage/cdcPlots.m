function f=cdcPlots(x)
%x=sortrows(x,[2,3]);

inds09=find(x(:,2)==2009);
inds10=find(x(:,2)==2010);
inds11=find(x(:,2)==2011);
x=x([inds09;inds10;inds11],:);
years=x(1,2):x(end,2);
%numyears=length(years);

twelve=1:12;
numyears=length(years);
y=[kron(years',ones(12,1)),repmat(twelve',numyears,1)];
y(:,3:17)=nan;
%Clunky way:
for i=1:numyears
    xi=x(x(:,2)==years(i),:);
    for j=1:12
        month=find(xi(:,3)==j);
        if isempty(month)==0
            yind=(i-1)*12+j;
            xij=nanmean(xi(month,:),1);
            y(yind,3:end)=xij(:,5:end);
        end
    end
end

incScaled=y(:,6)./y(:,7).*y(:,16);%.*y(:,end-1);
tvec1=y(:,1); tvec2=y(:,2);
tvec=1:size(incScaled,1);

tvec(isnan(incScaled))=[];
%tvec1(isnan(incScaled))=[];
tvec2(isnan(incScaled))=[];
incScaled(isnan(incScaled))=[];

jans=tvec(tvec2==1);%-.5;
ljans=length(jans);
fs=12; lw=2;
ticks=[1,13,25];
labs={'01/09','01/10','01/11'};
%
figure
hold on
mx=max(incScaled);
for j=1:ljans
    plot([jans(j),jans(j)],[0,mx],'k--','linewidth',lw)
end
plot(tvec,incScaled,'-','linewidth',lw)
%Manual (for plotting):
plot(tvec(5:12),incScaled(5:12),'-','linewidth',lw+1)
hold off
xlabel('Time (months)')
ylabel('Incidence')
xvec=1:6:tvec(end);
set(gca,'xtick',ticks,'xticklabel',labs)
set(gca,'fontsize',fs)
ax=gca;
ax.XAxis.MinorTickValues=tvec;
axis([0,tvec(end),0,mx])
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