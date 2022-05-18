%%
clc;
clear all;
%%
global N D LIMIT
N=50;       %colony size  
D=30;        %dimension
GENER=1000;
SAMP=10;
LIMIT=[-30,30];
PF=0.7;
SCALE1=100;
SCALE2=1;

HISTORY=zeros(GENER,SAMP);
%%  sample head
for sample=1:SAMP
%%  step 1
N_F=floor((0.9-rand*0.25)*N);       %number of female
N_M=N-N_F;      %number of male
%%  step 2
index_all=(1:1:N)';        %index of all spiders
index_female=(1:1:N_F)';       %index of female
index_male=(N_F+1:1:N)';       %index of male
mateRadius=(LIMIT(2)-LIMIT(1))/2./SCALE2;       %radius of mating
position_female=LIMIT(1)+rand(N_F,D)*(LIMIT(2)-LIMIT(1));
position_male=LIMIT(1)+rand(N_M,D)*(LIMIT(2)-LIMIT(1));
position_all=[position_female;position_male];
lastBest=inf;
%%  loop head
for generation=1:GENER
%%  step 3
fitness_female=f_Rosenbrock(position_female);
fitness_male=f_Rosenbrock(position_male);
fitness_all=[fitness_female;fitness_male];
fitness_globalBest=min(fitness_all);
fitness_globalWorst=max(fitness_all);
index_globalBest=find(fitness_all==fitness_globalBest,1);
index_globalWorst=find(fitness_all==fitness_globalWorst,1);
weight_female=(fitness_female-fitness_globalWorst)/...
    (fitness_globalBest-fitness_globalWorst);
weight_male=(fitness_male-fitness_globalWorst)/...
    (fitness_globalBest-fitness_globalWorst);
weight_all=[weight_female;weight_male];
if fitness_all(index_globalBest)<lastBest
    HISTORY(generation,sample)=fitness_all(index_globalBest);
else
    HISTORY(generation,sample)=lastBest;
end
lastBest=HISTORY(generation,sample);
%%   step 4
distance=zeros(N,N);
for i=1:N
    distance(:,i)=f_distance(position_all(i,:),position_all)./SCALE1;
    distance(i,i)=inf;
end
distance_copy=distance(1:1:N_F,:);

flag_higher=1;
count=1;
while sum(flag_higher)<N_F    
distance_closest=min(distance_copy,[],2);
index_closest=zeros(N_F,1);
for i=1:N_F
    index_closest(i)=find(distance_copy(i,:)==distance_closest(i),1);
end
flag_higher=weight_female<=weight_all(index_closest);
for i=1:N_F
    distance_copy(i,index_closest(i))=inf^(~flag_higher(i))*distance_copy(i,index_closest(i));
end
count=count+1;

if count==N_F+1 
%     index_flagHigher0=find(flag_higher==0); 
%     index_closest(index_flagHigher0)=index_flagHigher0;
%     distance_closest(index_flagHigher0)=0;
    flag_higher(flag_higher==0)=1;
end
end

Vibc=weight_all(index_closest).*exp(-distance_closest.^2);
Vibc(flag_higher==0)=0;
distance_globalBest=distance(:,index_globalBest);
Vibb=weight_all(index_globalBest).*exp(-distance_globalBest.^2);
Vibb(index_globalBest)=0;
R=rand(2);
pf=position_female;
position_closest=position_all(index_closest,:);
pc=position_closest(1:1:N_F,:);
position_globalBest=position_all(index_globalBest,:);
pb=position_globalBest;
vc=Vibc(1:1:N_F);
vb=Vibb(1:1:N_F);
flag_tend=rand(N_F,1)>PF;       % >PF minus sign & <PF plus sign
pf=pf+(-1).^flag_tend.*R(1).*vc.*(pc-pf)+(-1).^flag_tend.*R(2).*vb.*(pb-pf)+R(3)*(R(4)-0.5);
position_female=pf;
%%  step 5
weight_median=median(weight_male);
flag_aboveMedian=weight_male>=weight_median;

distance_male2female=distance(N_F+1:1:N,1:1:N_F);
distance_closestFemale=min(distance_male2female,[],2);
index_closestFemale=zeros(N_M,1);
for i=1:N_M
    index_closestFemale(i)=find(distance_male2female(i,:)==distance_closestFemale(i),1);
end
Vibf=weight_all(index_closestFemale).*exp(-distance_closestFemale.^2);

pm=position_male;
position_closestFemale=position_all(index_closestFemale,:);
pcf=position_closestFemale;
vf=Vibf;
SIGMA=sum(pm.*weight_male)./sum(weight_male);
pm=pm+flag_aboveMedian.*R(1).*vf.*(pcf-pm)+R(3)*(R(4)-0.5)...
    +(~flag_aboveMedian).*(SIGMA-pm);
position_male=pm;
%%  refresh fitness weight
% fitness_female=f_Rosenbrock(position_female);
% fitness_male=f_Rosenbrock(position_male);
% fitness_all=[fitness_female;fitness_male];
% fitness_globalBest=min(fitness_all);
% fitness_globalWorst=max(fitness_all);
% index_globalBest=find(fitness_all==fitness_globalBest,1);
% index_globalWorst=find(fitness_all==fitness_globalWorst,1);
% weight_female=(fitness_female-fitness_globalWorst)/...
%     (fitness_globalBest-fitness_globalWorst);
% weight_male=(fitness_male-fitness_globalWorst)/...
%     (fitness_globalBest-fitness_globalWorst);
% weight_all=[weight_female;weight_male];

%%  step 6
% weight_median=median(weight_male);
% flag_aboveMedian=weight_male>=weight_median;
index_dominant=index_male.*flag_aboveMedian;
index_dominant(index_dominant==0)=[];
position_dominant=position_all(index_dominant,:);
weight_dominant=weight_all(index_dominant);
distance_mate=distance(index_dominant,1:1:N_F);
flag_inMateRadius=distance_mate<mateRadius;

for i=1:size(index_dominant)
    num_mom=sum(flag_inMateRadius(i,:));
    if num_mom>0
        index_mom=(find(flag_inMateRadius(i,:)==1))';
        weight_mom=weight_female(index_mom);
        weight_mateGroup=[weight_mom;weight_dominant(i)];
        index_mateGroup=[index_mom;index_dominant(i)];
        p_chosen=weight_mateGroup./sum(weight_mateGroup);
        roulette=cumsum(p_chosen);
        dice=rand(1);
        index_chosen=find(dice<=roulette,1);
        index_newbrood=index_mateGroup(index_chosen);
        if weight_all(index_globalWorst)<weight_all(index_newbrood)
            position_all(index_globalWorst,:)=position_all(index_newbrood);
            weight_all(index_globalWorst)=weight_all(index_newbrood);
        end
    end
end

end
end
BEST_FITNESS=mean(HISTORY,2);

%%
figure
semilogy(mean(HISTORY,2)); %draw fitness 
figure
plot(mean(HISTORY,2));
tip=HISTORY(GENER,:);
AB=mean(tip);
MB=min(tip);
SD=std(tip);
DATA=[AB;MB;SD];





