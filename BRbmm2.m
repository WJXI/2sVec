function [Rbmm] = BRbmm2(Retraction,comp_state)
% comp_state=[01 02 12 012];
% Retraction is a bimodule map from bimodule A to B. comp_state is a bimodule C
% BRbmm2 is the relative bimodule map AC to BC, C in left
sigma=[1 0; 0 -1];
exchange=[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1];
[target_V, origin_V]=size(Retraction);
N_A=comp_state(1)+comp_state(2)+comp_state(3);
if N_A==3
    N_V=4;
elseif N_A==0
    N_V=1;
else
    N_V=2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N_V==4 % to modify
    Rbmm=zeros(target_V*2,origin_V*2);
    for i=1:origin_V
        caci=i-1;
        i_two=[fix(caci/4) fix((caci-4*fix(caci/4))/2) caci-4*fix(caci/4)-2*fix((caci-4*fix(caci/4))/2)];
        for j=1:target_V
            cacj=j-1;
            j_two=[fix(cacj/4) fix((cacj-4*fix(cacj/4))/2) cacj-4*fix(cacj/4)-2*fix((cacj-4*fix(cacj/4))/2)];
            for k=1:2
                new_i=[i_two(1:2) k-1 i_two(3)]*[8 4 2 1]';
                new_j=[j_two(1:2) k-1 j_two(3)]*[8 4 2 1]';
                Rbmm(new_j+1,new_i+1)=Retraction(j,i);
            end
        end
    end
    Rbmm=kron(eye(target_V/2),exchange)*kron(Retraction,eye(2))*kron(eye(origin_V/2),exchange);
elseif N_V==2 && comp_state(3)==0 % 3 config
    Rbmm=kron(Retraction,eye(2));
elseif N_V==2 && comp_state(3)==1 && N_A==1 % 1 config
    Rbmm=Retraction;
elseif N_V==2 && comp_state(3)==1 && comp_state(1)==1% 1 config
    Rbmm=Retraction;
elseif N_V==2 && comp_state(3)==1 && comp_state(2)==1  && comp_state(4)==0% 1 config with 1
    Rbmm=Retraction; 
elseif N_V==2 && comp_state(3)==1 && comp_state(2)==1  && comp_state(4)==1% 1 config with f
    Rbmm=kron(eye(target_V/2),sigma)*Retraction*kron(eye(origin_V/2),sigma); 
else
    Rbmm=Retraction;



end

