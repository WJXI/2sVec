function [BASISV,dimV] = Lbasis(object4,morphism4)
%LBASIS 
ob01=object4(1);
ob02=object4(2);
ob03=object4(3);
ob04=object4(4);
ob12=object4(5);
%ob13=object4(6);
ob14=object4(7);
ob23=object4(8);
ob24=object4(9);
ob34=object4(10);
mor012=morphism4(1);
%mor013=morphism4(2);
mor014=morphism4(3);
mor023=morphism4(4);
mor024=morphism4(5);
mor034=morphism4(6);
%mor123=morphism4(7);
mor124=morphism4(8);
%mor134=morphism4(9);
mor234=morphism4(10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Left path
[ARetraction1, ARetraction2]=lib_AssociBmm([ob02 ob03 ob04 ob23 ob24 ob34],[mor023 mor024 mor034 mor234]);
Acheck=sum(size(ARetraction2));
if Acheck==0
    zeta10=FRbmm(ARetraction1,[ob01 ob02 ob12 mor012]);
else
    zeta10=FRbmm(ARetraction1,[ob01 ob02 ob12 mor012]);
    zeta11=FRbmm(ARetraction2,[ob01 ob02 ob12 mor012]);
end
IC=interchange(object4,morphism4);
[BRetraction1, BRetraction2]=lib_AssociBmm([ob01 ob02 ob04 ob12 ob14 ob24],[mor012 mor014 mor024 mor124]);
Bcheck=sum(size(BRetraction2));
if Bcheck==0
    zeta20=FRbmm(BRetraction1,[ob23 ob24 ob34 mor234]);
else
    zeta20=FRbmm(BRetraction1,[ob23 ob24 ob34 mor234]);
    zeta21=FRbmm(BRetraction2,[ob23 ob24 ob34 mor234]);
end

if Acheck==0 && Bcheck==0
    Z00=zeta20*IC*zeta10;
    BASISV(1,:)=Z00(:);
    dimV=[1 1];
elseif Acheck==0 && Bcheck>0.1
    Z00=zeta20*IC*zeta10;
    Z01=zeta21*IC*zeta10;
    BASISV(1,:)=Z00(:);
    BASISV(2,:)=Z01(:);
    dimV=[1 2];
elseif Acheck>0.1 && Bcheck==0
    Z00=zeta20*IC*zeta10;
    Z10=zeta20*IC*zeta11;
    BASISV(1,:)=Z00(:);
    BASISV(2,:)=Z10(:);
    dimV=[2 1];
else
    Z00=zeta20*IC*zeta10;
    Z01=zeta21*IC*zeta10;
    Z10=zeta20*IC*zeta11;
    Z11=zeta21*IC*zeta11;
    BASISV(1,:)=Z00(:);
    BASISV(2,:)=Z01(:);
    BASISV(3,:)=Z10(:);
    BASISV(4,:)=Z11(:);
    dimV=[2 2];
end
end

