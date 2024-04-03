function [tj,dimVL,dimVR] = tj(object4,morphism4)
ob01=object4(1);
ob02=object4(2);
ob03=object4(3);
ob04=object4(4);
ob12=object4(5);
% ob13 J
ob14=object4(7);
ob23=object4(8);
ob24=object4(9);
ob34=object4(10);
mor012=morphism4(1); %P3
if mod(ob01+ob02+ob12,2)==0
    P3=mor012;
else
    P3=2;
end
% mor013 X
mor014=morphism4(3); %Q1
if mod(ob01+ob04+ob14,2)==0
    Q1=mor014;
else
    Q1=2;
end
mor023=morphism4(4); %P2
if mod(ob02+ob03+ob23,2)==0
    P2=mor023;
else
    P2=2;
end
% mor024 Z
mor034=morphism4(6); %P1
if mod(ob03+ob04+ob34,2)==0
    P1=mor034;
else
    P1=2;
end
% mor123 W
mor124=morphism4(8); %Q2
if mod(ob12+ob14+ob24,2)==0
    Q2=mor124;
else
    Q2=2;
end
% mor134 Y;
mor234=morphism4(10); %Q3
if mod(ob23+ob24+ob34,2)==0
    Q3=mor234;
else
    Q3=2;
end
[~,~,dim_Z,~] = three2two([P1 P2 P3],[Q1 Q2 Q3]);

%ifstate_object=[01 02 03 04 12 $13(J)$ 14 23 24 34]; $13(J)$ means 13(J) is undetermined 
%ifstate_morphism=[012 $013(X)$ 014 023 $024(Z)$ 034 $123(W)$ 124 $134(Y)$ 234]

if dim_Z==2
    [BASISV1,dimVL] = Lbasis(object4,morphism4);
    cac_morphism4=morphism4;
    cac_morphism4(5)=mod(morphism4(5)+1,2);
    [BASISV2] = Lbasis(object4,cac_morphism4);
    BASISV=[BASISV1;BASISV2];
else
    [BASISV,dimVL] = Lbasis(object4,morphism4);
end

[cac_Rstate,dimVR] = Rstate(object4,morphism4);


[Numbasis,~]=size(BASISV);
[Numstate,~]=size(cac_Rstate);
for i=1:Numbasis
    for j=1:Numstate
        tj(i,j)=cac_Rstate(j,:)*BASISV(i,:)';
    end
end
for i=1:size(tj,2)
    NORM=tj(:,i)'*tj(:,i);
    if NORM>0.01
        tj(:,i)=tj(:,i)/sqrt(NORM);
    end
end
if dim_Z==1
    tj=vpa(tj,2);
else
    tj=tj(1:1:dimVL(1)*dimVL(2),:);
    tj=vpa(tj,2);
end

end

