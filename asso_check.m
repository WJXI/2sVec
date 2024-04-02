function [output] = asso_check(state)
%check whether the configuration of four 1-morphisms in one associator satisfies fusion rule.
num_m=size(find(state==2));
if num_m(2)==1
    output=0;
elseif num_m(2)==3
    output=0;
elseif num_m(2)==0 && mod(state(1)+state(2)-state(3)-state(4),2)==1
    output=0;
else
    output=1;
end



end

