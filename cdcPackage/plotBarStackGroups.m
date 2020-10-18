function [] = plotBarStackGroups(stackData, groupLabels)
%% Plot a set of stacked bars, but group them according to labels provided.
%%
%% Params: 
%%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)
%%
%% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu
%%
%% 
NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);

cmap=lines(5);

% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
figure
    hold on; 
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    
    h(i,:) = bar(Y, 'stacked');
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
    
    for j=1:size(stackData,3)
        set(h(i,j),'FaceColor',cmap(j,:));
    end
end
hold off;
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);

%DH:
fs=10;
legString={'0-4','5-17','18-49','50-64','65+'};
ylabel('Hosp./1000')
set(gca,'fontsize',fs)
legend(legString,'location','NW')
box on
grid on
grid minor
end 
