function [Z,neworder_YWX,dim_Z,dim_YWX] = three2two(P,Q)
%initial state P1 P2 P3
%final state Q1 Q2 Q3
%path x is determined by Z
%path y is determined by X and W Y j.

P1=P(1);
P2=P(2);
P3=P(3);
Q1=Q(1);
Q2=Q(2);
Q3=Q(3);

% for path x
% P1 * P2 = Z * Q3      alpha
% Q1 * Q2 = Z * P3      beta
dim_Z=0;
for i=1:3 % i-1 represent Z
    state_1=[P1 P2 i-1 Q3]; % for alpha
    check_1=asso_check(state_1);
    state_2=[Q1 Q2 i-1 P3]; % for alpha
    check_2=asso_check(state_2);
    if sum(check_1+check_2)==2
        dim_Z=dim_Z+1;
        Z(dim_Z)=i-1;
    end
end

% for path y
% P2 * P3 = X * W       gamma
% P1 * X = Y * Q1       eta
% Y * W = Q2 * Q3       kappa
dim_YWX=0;
for k=1:3       %Y
    for j=1:3       %W
        for i=1:3       %X
            state_3=[P2 P3 i-1 j-1]; % for gamma
            state_4=[P1 i-1 k-1 Q1]; % for eta
            state_5=[j-1 k-1 Q2 Q3]; % for kappa
            check_3=asso_check(state_3);
            check_4=asso_check(state_4);
            check_5=asso_check(state_5);
            if sum(check_3+check_4+check_5)==3
                dim_YWX=dim_YWX+1;
                YWX(dim_YWX,:)=[k-1 j-1 i-1]; %Y W X
            end
        end
    end
end
for i=1:dim_YWX
    find_maj=find(YWX(i,:)>1);
    weight(i)=20*size(find_maj,2);
    cache_YWX=YWX(i,:);
    cache_YWX(find_maj)=[];
    if size(cache_YWX,2)==3
        weight(i)=weight(i)+4*cache_YWX(1)+2*cache_YWX(2)-cache_YWX(3);
    elseif size(cache_YWX,2)==2
        weight(i)=weight(i)+2*cache_YWX(1)-cache_YWX(2);
    elseif size(cache_YWX,2)==1
        weight(i)=weight(i)-cache_YWX(1);
    else
        weight(i)=weight(i);
    end
end
[~,order]=sort(weight);

neworder_YWX=zeros(dim_YWX,3);     % GIVE an order of the configuration YWX such that it is consistent with our article 
for i=1:dim_YWX
    neworder_YWX(dim_YWX+1-i,:)=YWX(order(i),:);
end



    
end

