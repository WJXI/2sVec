function [Rbmm] = BRbmm1(Retraction,comp_state)
% comp_state=[01 02 12 012];
% Retraction is a bimodule map from bimodule A to B. comp_state is a bimodule C
% BRbmm1 is the relative bimodule map AC to BC, C in right
sigma=[1 0; 0 -1];
AtoV=[1 0;0 -1i];
%AtoVf=[0 1i;1 0];
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
if N_V==4
    Rbmm=kron(Retraction,eye(2));
elseif N_V==2 && comp_state(1)==0 % 3 config
    Rbmm=kron(Retraction,eye(2));
elseif N_V==2 && comp_state(1)==1 && N_A==1 % 1 config
    Rbmm=Retraction;
elseif N_V==2 && comp_state(1)==1 && comp_state(3)==1% 1 config
    Rbmm=kron(eye(target_V/2),AtoV)*Retraction*kron(eye(origin_V/2),AtoV');
%elseif N_V==2 && comp_state(1)==1 && comp_state(3)==1 && comp_state(4)==1% 1 config with f
%    Rbmm=kron(eye(target_V/2),AtoV)*Retraction*kron(eye(origin_V/2),AtoV');
elseif N_V==2 && comp_state(1)==1 && comp_state(2)==1  && comp_state(4)==0% 1 config with 1
    Rbmm=Retraction; 
elseif N_V==2 && comp_state(1)==1 && comp_state(2)==1  && comp_state(4)==1% 1 config with f
    Rbmm=kron(eye(target_V/2),sigma)*Retraction*kron(eye(origin_V/2),sigma); 
else
    Rbmm=Retraction;



end

