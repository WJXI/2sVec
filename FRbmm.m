function [Rbmm] = FRbmm(Retraction,comp_state)
% comp_state=[01 02 12 012];
% Retraction is a bimodule map from bimodule A to B. comp_state is a bimodule C
% FRbmm is the relative bimodule map CA to CB

%[target_V, origin_V]=size(Retraction);
N_A=comp_state(1)+comp_state(2)+comp_state(3);
if N_A==3
    N_V=4;
elseif N_A==0
    N_V=1;
else
    N_V=2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N_V==4 
    Rbmm=kron(eye(2),Retraction);
elseif N_V==2 && comp_state(2)==0 
    Rbmm=kron(eye(2),Retraction);
elseif N_V==2 && comp_state(2)==1 
    Rbmm=Retraction;
else
    Rbmm=Retraction;



end

