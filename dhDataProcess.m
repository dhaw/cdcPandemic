function f=dhDataProcess(tab)

DV  = datevec(tab.OnsetDate);  % [N x 6] array
DV  = DV(:, 1:3);   % [N x 3] array, no time
DV2 = DV;
DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
%Result = cat(2, DV(:, 1), datenum(DV) - datenum(DV2));
%doy=result(:,2);
doyOnset = cat(1, datenum(DV) - datenum(DV2));

DV  = datevec(tab.AdmitDate);  % [N x 6] array
DV  = DV(:, 1:3);   % [N x 3] array, no time
DV2 = DV;
DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
%Result = cat(2, DV(:, 1), datenum(DV) - datenum(DV2));
%doy=result(:,2);
doyAdmit = cat(1, datenum(DV) - datenum(DV2));

doyOut=doyAdmit+tab.DaysHosp;

days=32:max([doyAdmit;doyOut]);
ldays=length(days);
inc=zeros(ldays,1);
prev=inc;
onset=inc;
doyOut(isnan(doyOut)==1)=days(end)+1;
for i=1:length(doyOnset)
    vi=doyAdmit(i)-31:doyOut(i)-32;
    inc(doyAdmit(i)-31)=inc(doyAdmit(i)-31)+1;
    prev(vi)=prev(vi)+1;
    onset(doyOnset(i)-31)=onset(doyOnset(i)-31)+1;
end

f=prev;

figure
cmap=lines(7);
fs=10; lw=2;
hold on
plot(days,cumsum(onset),'-','linewidth',lw,'color',cmap(2,:))
plot(days,cumsum(inc),'-','linewidth',lw,'color',cmap(1,:))
plot(days,inc,'--','linewidth',lw,'color',cmap(1,:))
plot(days,prev,'k-','linewidth',lw)
xlabel('Day (since 1st Feb)')
ylabel('Hospitalisations')
set(gca,'fontsize',fs)
legend('Total symptomatic','Total hospitalised','New hospitalisations','Number currently hospitalised','location','NW')
grid on
grid minor
box on