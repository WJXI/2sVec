clear all
clc
load('config18.mat'); %024=6 135=15
load('config024.mat');
load('config135.mat');
load('basis5_1.mat');
load('count.mat')
dim=count;
high_dim=find(dim>1);
[num_high_dim,~]=size(high_dim);
%check=zeros(1,num_high_dim);
for i=60019 %25 17
    state=high_dim(i);
    if mod(i,100)==0
        i/100
    end
    cac_024=unique(config_024(state,:));
    cac_024(cac_024>2)=[];
    cac_135=unique(config_135(state,:));
    cac_135(cac_135>2)=[];
    Map024=0;
    Map135=0;
    for j=1:size(cac_024,2)
        morphism5=[config(state,1:5) cac_024(j) config(state,6:13) cac_135(1) config(state,14:18)];
        object5=possibleobj(morphism5);
        morphism5=mod(morphism5,2);
        Map024=Map024+hexagon024(object5,morphism5);
    end
    for j=1:size(cac_135,2)
        morphism5=[config(state,1:5) cac_024(1) config(state,6:13) cac_135(j) config(state,14:18)];
        object5=possibleobj(morphism5);
        morphism5=mod(morphism5,2);
        Map135=Map135+hexagon135(object5,morphism5);
    end
    amp_find=find(abs(Map024)>0.2);
    amp(i)=Map024(amp_find(1))/Map135(amp_find(1));
    Map135=Map135*amp(i);
    check(i)=sum(sum(abs(Map024-Map135)));
end



