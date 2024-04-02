function [tjmatrix] = tjmatrix(ifstate_object, ifstate_morphism)
%ifstate is the configuration of initial and final state
%ifstate_object=[01 02 03 04 12 $13(J)$ 14 23 24 34]; $13(J)$ means 13(J) is undetermined 
%ifstate_morphism=[012 $013(X)$ 014 023 $024(Z)$ 034 $123(W)$ 124 $134(Y)$ 234]
%ifstate_object=[0 1 1 0 0 1 1 0 0]; %A M2 M1 K B N1 C N2 D[01 02 03 04 12 14 23 24 34]
%ifstate_morphism=[0 0 0 0 0 0]; %P3 Q1 P2 P1 Q2 Q3 [012 014 023 034 124 234]
ob01=ifstate_object(1);
ob02=ifstate_object(2);
ob03=ifstate_object(3);
ob04=ifstate_object(4);
ob12=ifstate_object(5);
% ob13 J
ob14=ifstate_object(6);
ob23=ifstate_object(7);
ob24=ifstate_object(8);
ob34=ifstate_object(9);
mor012=ifstate_morphism(1); %P3
if mod(ob01+ob02+ob12,2)==0
    P3=mor012;
else
    P3=2;
end
% mor013 X
mor014=ifstate_morphism(2); %Q1
if mod(ob01+ob04+ob14,2)==0
    Q1=mor014;
else
    Q1=2;
end
mor023=ifstate_morphism(3); %P2
if mod(ob02+ob03+ob23,2)==0
    P2=mor023;
else
    P2=2;
end
% mor024 Z
mor034=ifstate_morphism(4); %P1
if mod(ob03+ob04+ob34,2)==0
    P1=mor034;
else
    P1=2;
end
% mor123 W
mor124=ifstate_morphism(5); %Q2
if mod(ob12+ob14+ob24,2)==0
    Q2=mor124;
else
    Q2=2;
end
% mor134 Y;
mor234=ifstate_morphism(6); %Q3
if mod(ob23+ob24+ob34,2)==0
    Q3=mor234;
else
    Q3=2;
end
[Z,YWX,dim_Z,dim_YWX] = three2two([P1 P2 P3],[Q1 Q2 Q3]);

for i=1:dim_YWX
    morphism4(i,:)=[ifstate_morphism(1) mod(YWX(i,3),2) ifstate_morphism(2:3) mod(Z(1),2) ifstate_morphism(4)...
        mod(YWX(i,2),2) ifstate_morphism(5) mod(YWX(i,1),2) ifstate_morphism(6)];
    if YWX(i,3)==2
        ob13=mod(3-ob01-ob03,2);
    else
        ob13=mod(2-ob01-ob03,2);
    end
    object4(i,:)=[ifstate_object(1:5) ob13 ifstate_object(6:9)];
end

if dim_Z==2
    [BASISV1] = Lbasis(object4(1,:),morphism4(1,:));
    cac_morphism4=morphism4(1,:);
    cac_morphism4(5)=mod(Z(2),2);
    [BASISV2] = Lbasis(object4(1,:),cac_morphism4(1,:));
    BASISV=[BASISV1;BASISV2];
else
    [BASISV] = Lbasis(object4(1,:),morphism4(1,:));
end
totalRstate=[];
for i=1:dim_YWX
    [cac_Rstate] = Rstate(object4(i,:),morphism4(i,:));
    totalRstate=[totalRstate; cac_Rstate];
end

[Numbasis,~]=size(BASISV);
for i=1:Numbasis
    norm=BASISV(1,:)*BASISV(1,:)';
    if abs(norm)>0.01
        norm_basis(i,:)=BASISV(i,:)/sqrt(norm);
    else
        norm_basis(i,:)=BASISV(i,:);
    end
end

[Numstate,~]=size(totalRstate);
for i=1:Numstate
    %norm=totalRstate(i,:)*totalRstate(i,:)';
    if abs(norm)>0.01 
        norm_state(i,:)=totalRstate(i,:)/sqrt(norm);
    else
        norm_state(i,:)=totalRstate(i,:);
    end
end

for i=1:Numbasis
    for j=1:Numstate
        tj(i,j)=norm_state(j,:)*norm_basis(i,:)';
    end
end

tjmatrix=vpa(tj,2);
end

