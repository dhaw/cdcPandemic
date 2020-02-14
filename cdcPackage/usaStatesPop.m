function f=usaStatesPop(table)
these={'STATE','AGE','POPEST2010_CIV','WEEKLYRATE'};
table=table(:,ismember(table.Properties.VariableNames,these));

c=cell(10,1);%FIPS codes for states in each region
c{1}=[23,25,33,44,09];
c{2}=[36,34];
c{3}=[42,51,54,10,24,11];
c{4}=[21,37,45,47,28,01,13,12];
c{5}=[27,55,26,17,18,39];
c{6}=[35,48,40,04,22];
c{7}=[31,20,18,29];
c{8}=[30,56,49,08,38,46];
c{9}=[06,32,04,15];
c{10}=[53,16,41,02];

table.AGE(table.AGE>64)=-5;
table.AGE(table.AGE>24)=-4;
table.AGE(table.AGE>17)=-3;
table.AGE(table.AGE>4)=-2;
table.AGE(table.AGE>-1)=-1;
table.AGE=-table.AGE;

table=unstack(table,'POPEST2010_CIV','AGE');
table.REG=zeros(length(table.STATE),1);

for i=1:10
    ci=c{i};
    table.REG(ismember(table.STATE,ci)==1)=-i;
end
table(table.REG==0,:)=[];

table.REG=-table.REG;
table.STATE=[];
tab=table2array(table);
regs=tab(:,end);
tab(:,end)=[];
nbar=size(tab,2);
age=zeros(10,5);
for i=1:nbar
    age(:,i)=accumarray(regs,tab(:,i));
end
age=array2table(age);
age.Region=(1:10)';
age=age(:,[nbar+1,1:nbar]);
f=age;