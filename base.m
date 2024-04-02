object4=[0 1 1 0 0 1 1 1 0 0]; %  [01 02 03 04 12 13 14 23 24 34]
morphism4=[0 0 0 0 0 0 0 0 0 0];% [012 013 014 023 024 034 123 124 134 234]
ob01=object4(1);
ob02=object4(2);
ob03=object4(3);
ob04=object4(4);
ob12=object4(5);
ob13=object4(6);
ob14=object4(7);
ob23=object4(8);
ob24=object4(9);
ob34=object4(10);
mor012=morphism4(1);
mor013=morphism4(2);
mor014=morphism4(3);
mor023=morphism4(4);
mor024=morphism4(5);
mor034=morphism4(6);
mor123=morphism4(7);
mor124=morphism4(8);
mor134=morphism4(9);
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

[BRetraction1, BRetraction2]=lib_AssociBmm([ob01 ob02 ob04 ob12 ob14 ob24],[mor012 mor014 mor024 mor124]);
Bcheck=sum(size(BRetraction2));
if Bcheck==0
    zeta20=FRbmm(BRetraction1,[ob23 ob24 ob34 mor234]);
else
    zeta20=FRbmm(BRetraction1,[ob23 ob24 ob34 mor234]);
    zeta21=FRbmm(BRetraction2,[ob23 ob24 ob34 mor234]);
end

if Acheck==0 && Bcheck==0
    Z00=zeta20*zeta10;
    BASISV(1,:)=Z00(:);
    dimZ=1;
elseif Acheck==0 && Bcheck>0.1
    Z00=zeta20*zeta10;
    Z01=zeta21*zeta10;
    BASISV(1,:)=Z00(:);
    BASISV(2,:)=Z01(:);
    dimZ=2;
elseif Acheck>0.1 && Bcheck==0
    Z00=zeta20*zeta10;
    Z10=zeta20*zeta11;
    BASISV(1,:)=Z00(:);
    BASISV(2,:)=Z10(:);
    dimZ=2;
else
    Z00=zeta20*zeta10;
    Z01=zeta21*zeta10;
    Z10=zeta20*zeta11;
    Z11=zeta21*zeta11;
    BASISV(1,:)=Z00(:);
    BASISV(2,:)=Z01(:);
    BASISV(3,:)=Z10(:);
    BASISV(4,:)=Z11(:);
    dimZ=4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right path
%
%
%
%
%
%
[CRetraction1, CRetraction2]=lib_AssociBmm([ob01 ob02 ob03 ob12 ob13 ob23],[mor012 mor013 mor023 mor123]);
Ccheck=sum(size(CRetraction2));
if Ccheck==0
    zeta30=BRbmm1(CRetraction1,[ob03 ob04 ob34 mor034]);
else
    zeta30=BRbmm1(CRetraction1,[ob03 ob04 ob34 mor034]);
    zeta31=BRbmm1(CRetraction2,[ob03 ob04 ob34 mor034]);
end

[DRetraction1, DRetraction2]=lib_AssociBmm([ob01 ob03 ob04 ob13 ob14 ob34],[mor013,mor014,mor034,mor134]);
Dcheck=sum(size(DRetraction2));
if Dcheck==0
    zeta40=FRbmm(DRetraction1,[ob12 ob13 ob23 mor123]);
else
    zeta40=FRbmm(DRetraction1,[ob12 ob13 ob23 mor123]);
    zeta41=FRbmm(DRetraction2,[ob12 ob13 ob23 mor123]);
end

[ERetraction1, ERetraction2]=lib_AssociBmm([ob12 ob13 ob14 ob23 ob24 ob34],[mor123 mor124 mor134 mor234]);
Echeck=sum(size(ERetraction2));
if Echeck==0
    zeta50=BRbmm2(ERetraction1,[ob01 ob04 ob14 mor014]);
else
    zeta50=BRbmm2(ERetraction1,[ob01 ob04 ob14 mor014]);
    zeta51=BRbmm2(ERetraction2,[ob01 ob04 ob14 mor014]);
end

if Ccheck==0 && Dcheck==0 && Echeck==0
    YWX000=zeta50*zeta40*zeta30;
    Rstate(1,:)=YWX000(:);
    dimYWX=1;
elseif Ccheck==0 && Dcheck==0 && Echeck>0.1
    YWX000=zeta50*zeta40*zeta30;
    YWX001=zeta51*zeta40*zeta30;
    Rstate(1,:)=YWX000(:);
    Rstate(2,:)=YWX001(:);
    dimYWX=2;
elseif Ccheck==0 && Dcheck>0.1 && Echeck==0
    YWX000=zeta50*zeta40*zeta30;
    YWX010=zeta50*zeta41*zeta30;
    Rstate(1,:)=YWX000(:);
    Rstate(2,:)=YWX010(:);
    dimYWX=2;
elseif Ccheck>0.1 && Dcheck==0 && Echeck==0
    YWX000=zeta50*zeta40*zeta30;
    YWX100=zeta50*zeta40*zeta31;
    Rstate(1,:)=YWX000(:);
    Rstate(2,:)=YWX100(:);
    dimYWX=2;
elseif Ccheck==0 && Dcheck>0.1 && Echeck>0.1
    YWX000=zeta50*zeta40*zeta30;
    YWX001=zeta51*zeta40*zeta30;
    YWX010=zeta50*zeta41*zeta30;
    YWX011=zeta51*zeta41*zeta30;
    Rstate(1,:)=YWX000(:);
    Rstate(2,:)=YWX001(:);
    Rstate(3,:)=YWX010(:);
    Rstate(4,:)=YWX011(:);
    dimYWX=4;
elseif Ccheck>0.1 && Dcheck==0 && Echeck>0.1
    YWX000=zeta50*zeta40*zeta30;
    YWX001=zeta51*zeta40*zeta30;
    YWX100=zeta50*zeta40*zeta31;
    YWX101=zeta51*zeta40*zeta31;
    Rstate(1,:)=YWX000(:);
    Rstate(2,:)=YWX001(:);
    Rstate(3,:)=YWX100(:);
    Rstate(4,:)=YWX101(:);
    dimYWX=4;
elseif Ccheck>0.1 && Dcheck>0.1 && Echeck==0
    YWX000=zeta50*zeta40*zeta30;
    YWX010=zeta50*zeta41*zeta30;
    YWX100=zeta50*zeta40*zeta31;
    YWX110=zeta50*zeta41*zeta31;
    Rstate(1,:)=YWX000(:);
    Rstate(2,:)=YWX010(:);
    Rstate(3,:)=YWX100(:);
    Rstate(4,:)=YWX110(:);
    dimYWX=4;
else
    YWX000=zeta50*zeta40*zeta30;
    YWX001=zeta51*zeta40*zeta30;
    YWX010=zeta50*zeta41*zeta30;
    YWX011=zeta51*zeta41*zeta30;
    YWX100=zeta50*zeta40*zeta31;
    YWX101=zeta51*zeta40*zeta31;
    YWX110=zeta50*zeta41*zeta31;
    YWX111=zeta51*zeta41*zeta31;
    Rstate(1,:)=YWX000(:);
    Rstate(2,:)=YWX001(:);
    Rstate(3,:)=YWX010(:);
    Rstate(4,:)=YWX011(:);
    Rstate(5,:)=YWX100(:);
    Rstate(6,:)=YWX101(:);
    Rstate(7,:)=YWX110(:);
    Rstate(8,:)=YWX111(:);
    dimYWX=8;
end


