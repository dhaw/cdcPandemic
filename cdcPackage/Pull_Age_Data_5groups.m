% Code to get age-specific inputs for the modelling:
% - Hospitalisation data ('data') and populations ('popns')
% - Multipliers and ranges, accounting for coverage of flusurv-NET ('mults')
% - Mixing matrix, aggregated as required from Mossong data ('rho')

clear all;

% -------------------------------------------------------------------------
% --- FluSurv-NET data

[num, txt, raw] = xlsread('Age_data.xlsx','Hospital incidence');
mat = reshape(num,5,175/5,5);                                              % Dims: 1.Age group 2.Week 3.Field
popns = mat(:,1,4);
data = mat(:,:,3)';                                                        % Dims: 1.Week 2.Age group


% -------------------------------------------------------------------------
% --- Multipliers

[~,~,raw] = xlsread('Age_data.xlsx','Multipliers');
[a,b] = find(strcmp('Product',raw));
mults.mu = cell2mat(raw(a(2)+[2:6], b(2)))';                               % Mean for normal distribution on multipliers

% Aproximate sigma (standard deviation) using 2SD approximation
lohi = cell2mat(raw(a(2)+[2:6], [b(2)-1,b(2)+1]))';
mults.sig = diff(lohi,1)/4;                                                % Variance for normal distribution on multipliers


% -------------------------------------------------------------------------
% --- Set up the mixing matrix required for the age groups being modelled

% M_ij is the average daily number of contacts of type j --> i
[mixing, ~, ~] = xlsread('Age_data.xlsx','UK mixing matrix');

% Taken from Mossong data: note that of every three elements, first is wrt to census, second wrt sample, and third being odds ratio
tmp1 = [0.0568	0.0939	0.605	0.0597	0.1008	0.592	0.0641	0.1008	0.636	0.0656	0.1038	0.632	0.0642	0.0583	1.101	0.0616	0.0583	1.056	0.0703	0.0593	1.186	0.0777	0.0692	1.124	0.0759	0.0613	1.239	0.0665	0.0543	1.224	0.0614	0.0652	0.942	0.0648	0.0534	1.215	0.0511	0.0652	0.784	0.045	0.0267	1.688	0.1152	0.0296	3.887];
tmp2 = reshape(tmp1,3,length(tmp1)/3);
% Select numbers wrt census
age_props = tmp2(1,:);

% Find weights necessary for summing across columns

% % Age groups: 0-4, 5-24, 25-49, 50-64, 65+
% aggreg_inds = {1,  2:4,     5,  6:13, 14:15};  num_ages = length(aggreg_inds);

% Hosp Age groups: 0-4, 5-17, 18-49, 50-64, 65+
aggreg_inds =       {1,  2:4,  5:10, 11:13, 14:15};  num_ages = length(aggreg_inds);

% % ILI Age groups:  0-4, 5-24, 25-49, 50-64, 65+
% aggreg_inds =       {1,  2:5,  6:10, 11:13, 14:15};  num_ages = length(aggreg_inds);
 
age_weights = zeros(1,length(age_props));
for i = 1:length(aggreg_inds)
    inds = aggreg_inds{i};
    age_weights(inds) = age_props(inds)/sum(age_props(inds));
end
% Back to mixing matrix: aggregate downwards, just straight sums
interim = zeros(num_ages, size(mixing,1));
for i = 1:num_ages
   interim(i,:) = sum(mixing(aggreg_inds{i},:),1);
end
% Now do across columns, weighting by age
rho = zeros(num_ages);
for i = 1:num_ages
   rho(:,i) = sum(interim(:,aggreg_inds{i}).*repmat(age_weights(aggreg_inds{i}),size(interim,1),1),2);
end

save('data_5age_gps','data','popns','mults','rho'); 

return;

ff=figure; fs = 14; lw = 1.5;
plot(data,'linewidth',lw); 
xlabel('MMWR week','fontsize',fs); ylabel('EIP frequency','fontsize',fs);
xl = xlim; xl(1) = 0; xlim(xl);
legend('<5yo','5-17yo','18-24yo','25-64yo','>65yr');
set(ff,'Position',[440   373   396   425]);
