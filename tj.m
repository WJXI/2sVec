function [tj,dimVL,dimVR] = tj(object4,morphism4)
[BASISV,dimVL] = Lbasis(object4,morphism4);
[Rstate_cac,dimVR] = Rstate(object4,morphism4);

[Numbasis,~]=size(BASISV);
for i=1:Numbasis
    norm=BASISV(i,:)*BASISV(i,:)';
    if abs(norm)>0
        norm_basis(i,:)=BASISV(i,:)/sqrt(norm);
    else
        norm_basis(i,:)=BASISV(i,:);
    end
end

[Numstate,~]=size(Rstate_cac);
for i=1:Numstate
    norm=Rstate_cac(i,:)*Rstate_cac(i,:)';
    if abs(norm)>0
        norm_state(i,:)=Rstate_cac(i,:)/sqrt(norm);
    else
        norm_state(i,:)=Rstate_cac(i,:);
    end
end


for i=1:Numbasis
    for j=1:Numstate
        tj(i,j)=norm_state(j,:)*norm_basis(i,:)';
    end
end

tj=vpa(tj,2);
end

