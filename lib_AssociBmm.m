function [retraction1,retraction2] = lib_AssociBmm(object,morphism)
%object=[01 02 03 12 13 23] morphism=[012 013 023 123]
%This is the library of retractions.
%After fixing representative objects and 1-morphisms and quotient maps,
%here we list all used retractions of associator bimodule maps.
retraction2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1
AftoA=[0 1; 1 0];
fAtoA=[0 1; -1 0];
fAAf=fAtoA'*AftoA;
AtoV=[1 0; 0 -1i];
VtoVf=[0 -1; 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2
toA=[1 0 0 1;0 1 1 0];
tofA=[0 1 -1 0; 1 0 0 -1];
idA=[1 0 0 0; 0 1 0 0];
idfA=[0 0 1 0; 0 0 0 1];
toV1=[1 0 0 1i;0 1 -1i 0];
toVf=[0 1 1i 0; 1 0 0 -1i];
V1to=[1 0 0 -1i;0 -1i 1 0]';
Vfto=[0 1i 1 0;1 0 0 1i]';
id2=eye(2);
% merge
VftofV=[-1 0; 0 1];
V1tofVf=[1 0; 0 -1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3-0
AVtoAV30=[1 0 0 1i; 0 1 -1i 0; 0 1i -1 0; -1i 0 0 -1];
AVtoVA30=[1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 -1];
VAtoAV30=[1 0 0 -1; 0 1 1 0; 0 -1 1 0; 1 0 0 1];
% merge
AVftoAV30=[0 0 1 0; 0 0 0 -1; 1 0 0 0; 0 -1 0 0];
VfAtoVA30=[0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2-1
VAtoAA21=[1 0 0 1; 0 1 1 0; 0 -1i 1i 0; -1i 0 0 1i];
AAtoAAt21=[1 0 0 1; 0 1 1 0; 0 1 -1 0; 1 0 0 -1];
VAtoAAt21=AAtoAAt21*VAtoAA21/(1-1i);
id21=eye(4);
% merge
AfAtoAA21=[0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0];
fAAtoAA21=[0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0];
AfAttoAAt21=[0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0];
VfAtoVA21=[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4
VAtoVA31=[1 0 0 1i; 0 1 1i 0; 0 1 -1i 0; 1 0 0 -1i];
VfAtoVfA31=[1 0 0 -1i; 0 1 -1i 0; 0 -1 -1i 0; -1 0 0 -1i];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class AAA-A [1 / 1 1 / 1]
if sum(abs(object-[1 1 1 1 1 1]))==0 && sum(abs(morphism-[0 0 0 0]))==0
    retraction1 = kron(V1to*toV1,eye(2));
    retraction2 = kron(Vfto*toVf,eye(2));

elseif sum(abs(object-[1 1 1 1 0 1]))==0 && mod(sum(morphism),2)==0
    retraction1 = VAtoVA31*kron(toV1,eye(2));
elseif sum(abs(object-[1 1 1 1 0 1]))==0 && mod(sum(morphism),2)==1
    retraction1 = VfAtoVfA31*kron(toVf,eye(2));

elseif sum(abs(object-[1 0 1 1 1 1]))==0 && mod(sum(morphism),2)==0
    retraction1 = kron(V1to,eye(2));
elseif sum(abs(object-[1 0 1 1 1 1]))==0 && mod(sum(morphism),2)==1
    retraction1 = kron(Vfto,eye(2));

elseif sum(abs(object-[1 0 1 1 0 1]))==0 && mod(morphism(1)+morphism(3),2)==0
    retraction1 = VAtoVA31;
elseif sum(abs(object-[1 0 1 1 0 1]))==0 && mod(morphism(1)+morphism(3),2)==1
    retraction1 = VfAtoVfA31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class AAA-1 [1 / 0 1 / 1]
elseif sum(abs(object-[1 1 0 1 1 1]))==0 && morphism(2)+morphism(3)==0
    retraction1 = AVtoAV30;
elseif sum(abs(object-[1 1 0 1 1 1]))==0 && mod(morphism(2)+morphism(3),2)==1 && morphism(3)==1
    retraction1 = AVtoAV30*AVftoAV30;
elseif sum(abs(object-[1 1 0 1 1 1]))==0 && mod(morphism(2)+morphism(3),2)==1 && morphism(2)==1
    retraction1 = AVftoAV30'*AVtoAV30; 
elseif sum(abs(object-[1 1 0 1 1 1]))==0 && morphism(2)==1 && morphism(3)==1
    retraction1 = AVftoAV30'*AVtoAV30*AVftoAV30;

elseif sum(abs(object-[1 1 0 1 0 1]))==0 && morphism(3)==0 && morphism(4)==0
    retraction1=AVtoVA30;
elseif sum(abs(object-[1 1 0 1 0 1]))==0 && morphism(3)==0 && morphism(4)==1
    retraction1=VfAtoVA30'*AVtoVA30;
elseif sum(abs(object-[1 1 0 1 0 1]))==0 && morphism(3)==1 && morphism(4)==0
    retraction1=AVtoVA30*AVftoAV30;
elseif sum(abs(object-[1 1 0 1 0 1]))==0 && morphism(3)==1 && morphism(4)==1
    retraction1=VfAtoVA30'*AVtoVA30*AVftoAV30;

elseif sum(abs(object-[1 0 0 1 1 1]))==0 && morphism(1)==0 && morphism(2)==0
    retraction1=VAtoAV30;
elseif sum(abs(object-[1 0 0 1 1 1]))==0 && morphism(1)==0 && morphism(2)==1
    retraction1=AVftoAV30'*VAtoAV30;
elseif sum(abs(object-[1 0 0 1 1 1]))==0 && morphism(1)==1 && morphism(2)==0
    retraction1=VAtoAV30*VfAtoVA30;
elseif sum(abs(object-[1 0 0 1 1 1]))==0 && morphism(1)==1 && morphism(2)==1
    retraction1=AVftoAV30'*VAtoAV30*VfAtoVA30;

elseif sum(abs(object-[1 0 0 1 0 1]))==0 && morphism(1)==0 && morphism(4)==0
    retraction1=VAtoVA31;
elseif sum(abs(object-[1 0 0 1 0 1]))==0 && morphism(1)==0 && morphism(4)==1
    retraction1=VfAtoVA30'*VAtoVA31;
elseif sum(abs(object-[1 0 0 1 0 1]))==0 && morphism(1)==1 && morphism(4)==0
    retraction1=VAtoVA31*VfAtoVA30;
elseif sum(abs(object-[1 0 0 1 0 1]))==0 && morphism(1)==1 && morphism(4)==1
    retraction1=VfAtoVfA31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class AA1-A [1 / 1 1 / 0]
elseif sum(abs(object-[1 1 1 1 1 0]))==0 && morphism(3)+morphism(4)==0
    retraction1=id21;
elseif sum(abs(object-[1 1 1 1 1 0]))==0 && morphism(3)==1 && morphism(4)==0
    retraction1=id21*AfAtoAA21;
elseif sum(abs(object-[1 1 1 1 1 0]))==0 && morphism(3)==0 && morphism(4)==1
    retraction1=fAAtoAA21'*id21;
elseif sum(abs(object-[1 1 1 1 1 0]))==0 && morphism(3)==1 && morphism(4)==1
    retraction1=fAAtoAA21'*id21*AfAtoAA21;

elseif sum(abs(object-[1 1 1 1 0 0]))==0 && morphism(3)+morphism(2)==0
    retraction1=AAtoAAt21;
elseif sum(abs(object-[1 1 1 1 0 0]))==0 && morphism(3)==0 && morphism(2)==1
    retraction1=AfAttoAAt21'*AAtoAAt21;
elseif sum(abs(object-[1 1 1 1 0 0]))==0 && morphism(3)==1 && morphism(2)==0
    retraction1=AAtoAAt21*AfAtoAA21;
elseif sum(abs(object-[1 1 1 1 0 0]))==0 && morphism(3)==1 && morphism(2)==1
    retraction1=AfAttoAAt21'*AAtoAAt21*AfAtoAA21;

elseif sum(abs(object-[1 0 1 1 1 0]))==0 && morphism(1)+morphism(4)==0
    retraction1=VAtoAA21;
elseif sum(abs(object-[1 0 1 1 1 0]))==0 && morphism(1)==0 && morphism(4)==1
    retraction1=fAAtoAA21'*VAtoAA21;
elseif sum(abs(object-[1 0 1 1 1 0]))==0 && morphism(1)==1 && morphism(4)==0
    retraction1=VAtoAA21*VfAtoVA21;
elseif sum(abs(object-[1 0 1 1 1 0]))==0 && morphism(1)==1 && morphism(4)==1
    retraction1=fAAtoAA21'*VAtoAA21*VfAtoVA21;

elseif sum(abs(object-[1 0 1 1 0 0]))==0 && morphism(1)==0 && morphism(2)==0
    retraction1=VAtoAAt21;
elseif sum(abs(object-[1 0 1 1 0 0]))==0 && morphism(1)==0 && morphism(2)==1
    retraction1=AfAttoAAt21'*VAtoAAt21;
elseif sum(abs(object-[1 0 1 1 0 0]))==0 && morphism(1)==1 && morphism(2)==0
    retraction1=VAtoAAt21*VfAtoVA21;
elseif sum(abs(object-[1 0 1 1 0 0]))==0 && morphism(1)==1 && morphism(2)==1
    retraction1=AfAttoAAt21'*VAtoAAt21*VfAtoVA21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class A1A-A [1 / 1 0 / 1]
elseif sum(abs(object-[1 1 1 0 1 1]))==0 && morphism(1)==0 && morphism(4)==0
    retraction1=id21;
elseif sum(abs(object-[1 1 1 0 1 1]))==0 && morphism(1)==0 && morphism(4)==1
    retraction1=fAAtoAA21'*id21;
elseif sum(abs(object-[1 1 1 0 1 1]))==0 && morphism(1)==1 && morphism(4)==0
    retraction1=id21*fAAtoAA21;
elseif sum(abs(object-[1 1 1 0 1 1]))==0 && morphism(1)==1 && morphism(4)==1
    retraction1=id21;

elseif sum(abs(object-[1 1 1 0 0 1]))==0 && morphism(1)==0 && morphism(2)==0
    retraction1=AAtoAAt21;
elseif sum(abs(object-[1 1 1 0 0 1]))==0 && morphism(1)==0 && morphism(2)==1
    retraction1=AfAttoAAt21'*AAtoAAt21;
elseif sum(abs(object-[1 1 1 0 0 1]))==0 && morphism(1)==1 && morphism(2)==0
    retraction1=AAtoAAt21*fAAtoAA21;
elseif sum(abs(object-[1 1 1 0 0 1]))==0 && morphism(1)==1 && morphism(2)==1
    retraction1=AfAttoAAt21'*AAtoAAt21*fAAtoAA21;

elseif sum(abs(object-[1 0 1 0 1 1]))==0 && morphism(3)==0 && morphism(4)==0
    retraction1=id21;
elseif sum(abs(object-[1 0 1 0 1 1]))==0 && morphism(3)==0 && morphism(4)==1
    retraction1=fAAtoAA21'*id21;
elseif sum(abs(object-[1 0 1 0 1 1]))==0 && morphism(3)==1 && morphism(4)==0
    retraction1=id21*AfAtoAA21;
elseif sum(abs(object-[1 0 1 0 1 1]))==0 && morphism(3)==1 && morphism(4)==1
    retraction1=fAAtoAA21'*id21*AfAtoAA21;

elseif sum(abs(object-[1 0 1 0 0 1]))==0 && morphism(3)==0 && morphism(2)==0
    retraction1=AAtoAAt21;
elseif sum(abs(object-[1 0 1 0 0 1]))==0 && morphism(3)==0 && morphism(2)==1
    retraction1=AfAttoAAt21'*AAtoAAt21;
elseif sum(abs(object-[1 0 1 0 0 1]))==0 && morphism(3)==1 && morphism(2)==0
    retraction1=AAtoAAt21*AfAtoAA21;
elseif sum(abs(object-[1 0 1 0 0 1]))==0 && morphism(3)==1 && morphism(2)==1
    retraction1=AfAttoAAt21'*AAtoAAt21*AfAtoAA21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 1AA-A [0 / 1 1 / 1]
elseif sum(abs(object-[0 1 1 1 1 1]))==0 && morphism(1)==0 && morphism(2)==0
    retraction1=id21;
elseif sum(abs(object-[0 1 1 1 1 1]))==0 && morphism(1)==0 && morphism(2)==1
    retraction1=AfAtoAA21'*id21;
elseif sum(abs(object-[0 1 1 1 1 1]))==0 && morphism(1)==1 && morphism(2)==0
    retraction1=id21*fAAtoAA21;
elseif sum(abs(object-[0 1 1 1 1 1]))==0 && morphism(1)==1 && morphism(2)==1
    retraction1=AfAtoAA21'*id21*fAAtoAA21;

elseif sum(abs(object-[0 1 1 1 0 1]))==0 && morphism(1)==0 && morphism(4)==0
    retraction1=VAtoAA21';
elseif sum(abs(object-[0 1 1 1 0 1]))==0 && morphism(1)==0 && morphism(4)==1
    retraction1=VfAtoVA21'*VAtoAA21';
elseif sum(abs(object-[0 1 1 1 0 1]))==0 && morphism(1)==1 && morphism(4)==0
    retraction1=VAtoAA21'*fAAtoAA21;
elseif sum(abs(object-[0 1 1 1 0 1]))==0 && morphism(1)==1 && morphism(4)==1
    retraction1=VfAtoVA21'*VAtoAA21'*fAAtoAA21;

elseif sum(abs(object-[0 0 1 1 1 1]))==0 && morphism(3)==0 && morphism(2)==0
    retraction1=id21;
elseif sum(abs(object-[0 0 1 1 1 1]))==0 && morphism(3)==0 && morphism(2)==1
    retraction1=AfAtoAA21'*id21;
elseif sum(abs(object-[0 0 1 1 1 1]))==0 && morphism(3)==1 && morphism(2)==0
    retraction1=id21*AfAtoAA21;
elseif sum(abs(object-[0 0 1 1 1 1]))==0 && morphism(3)==1 && morphism(2)==1
    retraction1=AfAtoAA21'*id21*AfAtoAA21;

elseif sum(abs(object-[0 0 1 1 0 1]))==0 && morphism(3)==0 && morphism(4)==0
    retraction1=VAtoAA21';
elseif sum(abs(object-[0 0 1 1 0 1]))==0 && morphism(3)==0 && morphism(4)==1
    retraction1=VfAtoVA21'*VAtoAA21';
elseif sum(abs(object-[0 0 1 1 0 1]))==0 && morphism(3)==1 && morphism(4)==0
    retraction1=VAtoAA21'*AfAtoAA21;
elseif sum(abs(object-[0 0 1 1 0 1]))==0 && morphism(3)==1 && morphism(4)==1
    retraction1=VfAtoVA21'*VAtoAA21'*AfAtoAA21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class AA1-1 [1 / 0 1 / 0]
elseif sum(abs(object-[1 1 0 1 1 0]))==0 && morphism(2)==0 && morphism(4)==0
    retraction1=toV1;
elseif sum(abs(object-[1 1 0 1 1 0]))==0 && morphism(2)==0 && morphism(4)==1
    retraction1=VftofV*toVf;
elseif sum(abs(object-[1 1 0 1 1 0]))==0 && morphism(2)==1 && morphism(4)==0
    retraction1=toVf;
elseif sum(abs(object-[1 1 0 1 1 0]))==0 && morphism(2)==1 && morphism(4)==1
    retraction1=V1tofVf*toV1;

elseif sum(abs(object-[1 1 0 1 0 0]))==0
    retraction1=V1to*toV1;
    retraction2=Vfto*toVf;

elseif sum(abs(object-[1 0 0 1 1 0]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(4)==0
    retraction1=id2;
elseif sum(abs(object-[1 0 0 1 1 0]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(4)==1
    retraction1=V1tofVf;
elseif sum(abs(object-[1 0 0 1 1 0]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(4)==0
    retraction1=id2;
elseif sum(abs(object-[1 0 0 1 1 0]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(4)==1
    retraction1=VftofV;

elseif sum(abs(object-[1 0 0 1 0 0]))==0 && mod(morphism(1)+morphism(3),2)==0
    retraction1=V1to;
elseif sum(abs(object-[1 0 0 1 0 0]))==0 && mod(morphism(1)+morphism(3),2)==1
    retraction1=Vfto;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class A1A-1 [1 / 0 0 / 1]
elseif sum(abs(object-[1 1 0 0 1 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && mod(morphism(2)+morphism(3),2)==0
    retraction1=id2;
elseif sum(abs(object-[1 1 0 0 1 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && mod(morphism(2)+morphism(3),2)==1
    retraction1=V1tofVf;
elseif sum(abs(object-[1 1 0 0 1 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && mod(morphism(2)+morphism(3),2)==0
    retraction1=id2;
elseif sum(abs(object-[1 1 0 0 1 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && mod(morphism(2)+morphism(3),2)==1
    retraction1=VftofV;

elseif sum(abs(object-[1 1 0 0 0 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(3)==0
    retraction1=V1to;
elseif sum(abs(object-[1 1 0 0 0 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(3)==1
    retraction1=V1to*V1tofVf';
elseif sum(abs(object-[1 1 0 0 0 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(3)==0
    retraction1=Vfto*VftofV';
elseif sum(abs(object-[1 1 0 0 0 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(3)==1
    retraction1=Vfto;

elseif sum(abs(object-[1 0 0 0 1 1]))==0 && mod(morphism(2)+morphism(4),2)==0 && morphism(2)==0
    retraction1=toV1;
elseif sum(abs(object-[1 0 0 0 1 1]))==0 && mod(morphism(2)+morphism(4),2)==0 && morphism(2)==1
    retraction1=V1tofVf*toV1;
elseif sum(abs(object-[1 0 0 0 1 1]))==0 && mod(morphism(2)+morphism(4),2)==1 && morphism(2)==0
    retraction1=VftofV*toVf;
elseif sum(abs(object-[1 0 0 0 1 1]))==0 && mod(morphism(2)+morphism(4),2)==1 && morphism(2)==1
    retraction1=toVf;

elseif sum(abs(object-[1 0 0 0 0 1]))==0
    retraction1=V1to*toV1;
    retraction2=Vfto*toVf;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 1AA-1 [0 / 0 1 / 1]
elseif sum(abs(object-[0 1 0 1 1 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(3)==0
    retraction1=toV1';
elseif sum(abs(object-[0 1 0 1 1 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(3)==1
    retraction1=toV1'*V1tofVf';
elseif sum(abs(object-[0 1 0 1 1 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(3)==0
    retraction1=toVf'*VftofV';
elseif sum(abs(object-[0 1 0 1 1 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(3)==1
    retraction1=toVf';

elseif sum(abs(object-[0 1 0 1 0 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(3)==0
    retraction1=id2;
elseif sum(abs(object-[0 1 0 1 0 1]))==0 && mod(morphism(1)+morphism(3),2)==0 && morphism(3)==1
    retraction1=id2*V1tofVf';
elseif sum(abs(object-[0 1 0 1 0 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(3)==0
    retraction1=id2*VftofV';
elseif sum(abs(object-[0 1 0 1 0 1]))==0 && mod(morphism(1)+morphism(3),2)==1 && morphism(3)==1
    retraction1=id2;

elseif sum(abs(object-[0 0 0 1 1 1]))==0 
    retraction1=toV1'*toV1;
    retraction2=toVf'*toVf;

elseif sum(abs(object-[0 0 0 1 0 1]))==0 && mod(morphism(2)+morphism(4),2)==0
    retraction1=toV1;
elseif sum(abs(object-[0 0 0 1 0 1]))==0 && mod(morphism(2)+morphism(4),2)==1
    retraction1=toVf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 11A-A [0 / 1 0 / 1]
elseif sum(abs(object-[0 1 1 0 1 1]))==0 && mod(morphism(2)+morphism(4),2)==0
    retraction1=idA;
elseif sum(abs(object-[0 1 1 0 1 1]))==0 && mod(morphism(2)+morphism(4),2)==1
    retraction1=idfA;

elseif sum(abs(object-[0 1 1 0 0 1]))==0 
    retraction1=toA'*idA;
    retraction2=tofA'*idfA;

elseif sum(abs(object-[0 0 1 0 1 1]))==0
    retraction1=id2;
    
elseif sum(abs(object-[0 0 1 0 0 1]))==0 && mod(morphism(1)+morphism(3),2)==0
    retraction1=toA';
elseif sum(abs(object-[0 0 1 0 0 1]))==0 && mod(morphism(1)+morphism(3),2)==1
    retraction1=tofA';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 1A1-A [0 / 1 1 / 0]
elseif sum(abs(object-[0 1 1 1 1 0]))==0
    retraction1=id2;

elseif sum(abs(object-[0 1 1 1 0 0]))==0 && mod(morphism(1)+morphism(3),2)==0
    retraction1=toA';
elseif sum(abs(object-[0 1 1 1 0 0]))==0 && mod(morphism(1)+morphism(3),2)==1
    retraction1=tofA';

elseif sum(abs(object-[0 0 1 1 1 0]))==0 && mod(morphism(2)+morphism(4),2)==0
    retraction1=toA;
elseif sum(abs(object-[0 0 1 1 1 0]))==0 && mod(morphism(2)+morphism(4),2)==1
    retraction1=tofA;

elseif sum(abs(object-[0 0 1 1 0 0]))==0 
    retraction1=toA'*toA;
    retraction2=tofA'*tofA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class A11-A [1 / 1 0 / 0]
elseif sum(abs(object-[1 1 1 0 1 0]))==0 && mod(morphism(1)+morphism(3),2)==0
    retraction1=toA';
elseif sum(abs(object-[1 1 1 0 1 0]))==0 && mod(morphism(1)+morphism(3),2)==1
    retraction1=tofA';

elseif sum(abs(object-[1 1 1 0 0 0]))==0
    retraction1=id2;

elseif sum(abs(object-[1 0 1 0 1 0]))==0
    retraction1=toA'*toA;
    retraction2=tofA'*tofA;

elseif sum(abs(object-[1 0 1 0 0 0]))==0 && mod(morphism(2)+morphism(4),2)==0
    retraction1=toA;
elseif sum(abs(object-[1 0 1 0 0 0]))==0 && mod(morphism(2)+morphism(4),2)==1
    retraction1=tofA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 111-A [0 / 1 0 / 0]
elseif sum(abs(object-[0 1 1 0 1 0]))==0 && mod(morphism(2)+morphism(3),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 1 1 0 1 0]))==0 && mod(morphism(2)+morphism(3),2)==1
    retraction1=AftoA;

elseif sum(abs(object-[0 1 1 0 0 0]))==0 && mod(morphism(4)+morphism(3),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 1 1 0 0 0]))==0 && mod(morphism(4)+morphism(3),2)==1
    retraction1=AftoA;

elseif sum(abs(object-[0 0 1 0 1 0]))==0 && mod(morphism(2)+morphism(1),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 0 1 0 1 0]))==0 && mod(morphism(2)+morphism(1),2)==1
    retraction1=AftoA;

elseif sum(abs(object-[0 0 1 0 0 0]))==0 && mod(morphism(4)+morphism(1),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 0 1 0 0 0]))==0 && mod(morphism(4)+morphism(1),2)==1
    retraction1=AftoA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 11A-1 [0 / 0 0 / 1]
elseif sum(abs(object-[0 1 0 0 1 1]))==0 && morphism(3)==0 && morphism(4)==0
    retraction1=id2;
elseif sum(abs(object-[0 1 0 0 1 1]))==0 && morphism(3)==0 && morphism(4)==1
    retraction1=fAtoA';
elseif sum(abs(object-[0 1 0 0 1 1]))==0 && morphism(3)==1 && morphism(4)==0
    retraction1=AftoA;
elseif sum(abs(object-[0 1 0 0 1 1]))==0 && morphism(3)==1 && morphism(4)==1
    retraction1=fAAf;

elseif sum(abs(object-[0 1 0 0 0 1]))==0 && mod(morphism(2)+morphism(3),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 1 0 0 0 1]))==0 && mod(morphism(2)+morphism(3),2)==1
    retraction1=AftoA;

elseif sum(abs(object-[0 0 0 0 1 1]))==0 && mod(morphism(1)+morphism(4),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 0 0 0 1 1]))==0 && morphism(1)==1 && morphism(4)==0
    retraction1=fAtoA;
elseif sum(abs(object-[0 0 0 0 1 1]))==0 && morphism(1)==0 && morphism(4)==1
    retraction1=fAtoA';

elseif sum(abs(object-[0 0 0 0 0 1]))==0 && morphism(1)==0 && morphism(2)==0
    retraction1=id2;
elseif sum(abs(object-[0 0 0 0 0 1]))==0 && morphism(1)==0 && morphism(2)==1
    retraction1=AftoA;
elseif sum(abs(object-[0 0 0 0 0 1]))==0 && morphism(1)==1 && morphism(2)==0
    retraction1=fAtoA;
elseif sum(abs(object-[0 0 0 0 0 1]))==0 && morphism(1)==1 && morphism(2)==1
    retraction1=fAAf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 1A1-1 [0 / 0 1 / 0]
elseif sum(abs(object-[0 1 0 1 1 0]))==0 && mod(morphism(1)+morphism(4),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 1 0 1 1 0]))==0 && morphism(1)==1 && morphism(4)==0
    retraction1=fAtoA;
elseif sum(abs(object-[0 1 0 1 1 0]))==0 && morphism(1)==0 && morphism(4)==1
    retraction1=fAtoA';

elseif sum(abs(object-[0 1 0 1 0 0]))==0 && morphism(1)==0 && morphism(2)==0
    retraction1=id2;
elseif sum(abs(object-[0 1 0 1 0 0]))==0 && morphism(1)==0 && morphism(2)==1
    retraction1=AftoA;
elseif sum(abs(object-[0 1 0 1 0 0]))==0 && morphism(1)==1 && morphism(2)==0
    retraction1=fAtoA;
elseif sum(abs(object-[0 1 0 1 0 0]))==0 && morphism(1)==1 && morphism(2)==1
    retraction1=fAAf;

elseif sum(abs(object-[0 0 0 1 1 0]))==0 && morphism(3)==0 && morphism(4)==0
    retraction1=id2;
elseif sum(abs(object-[0 0 0 1 1 0]))==0 && morphism(3)==0 && morphism(4)==1
    retraction1=fAtoA';
elseif sum(abs(object-[0 0 0 1 1 0]))==0 && morphism(3)==1 && morphism(4)==0
    retraction1=AftoA;
elseif sum(abs(object-[0 0 0 1 1 0]))==0 && morphism(3)==1 && morphism(4)==1
    retraction1=fAAf;

elseif sum(abs(object-[0 0 0 1 0 0]))==0 && mod(morphism(3)+morphism(2),2)==0
    retraction1=id2;
elseif sum(abs(object-[0 0 0 1 0 0]))==0 && mod(morphism(3)+morphism(2),2)==1
    retraction1=AftoA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class A11-1 [1 / 0 0 / 0]
elseif sum(abs(object-[1 1 0 0 1 0]))==0 && morphism(1)==0 && morphism(2)==0
    retraction1=AtoV;
elseif sum(abs(object-[1 1 0 0 1 0]))==0 && morphism(1)==0 && morphism(2)==1
    retraction1=VtoVf*AtoV;
elseif sum(abs(object-[1 1 0 0 1 0]))==0 && morphism(1)==1 && morphism(2)==0
    retraction1=AtoV*fAtoA;
elseif sum(abs(object-[1 1 0 0 1 0]))==0 && morphism(1)==1 && morphism(2)==1
    retraction1=VtoVf*AtoV*fAtoA;

elseif sum(abs(object-[1 1 0 0 0 0]))==0 && mod(morphism(1)+morphism(4),2)==0
    retraction1=id2;
elseif sum(abs(object-[1 1 0 0 0 0]))==0 && morphism(1)==1 && morphism(4)==0
    retraction1=fAtoA;
elseif sum(abs(object-[1 1 0 0 0 0]))==0 && morphism(1)==0 && morphism(4)==1
    retraction1=fAtoA';
    
elseif sum(abs(object-[1 0 0 0 1 0]))==0 && morphism(3)==0 && morphism(2)==0
    retraction1=AtoV;
elseif sum(abs(object-[1 0 0 0 1 0]))==0 && morphism(3)==0 && morphism(2)==1
    retraction1=VtoVf*AtoV;
elseif sum(abs(object-[1 0 0 0 1 0]))==0 && morphism(3)==1 && morphism(2)==0
    retraction1=AtoV*AftoA;
elseif sum(abs(object-[1 0 0 0 1 0]))==0 && morphism(3)==1 && morphism(2)==1
    retraction1=VtoVf*AtoV*AftoA;

elseif sum(abs(object-[1 0 0 0 0 0]))==0 && morphism(3)==0 && morphism(4)==0
    retraction1=id2;
elseif sum(abs(object-[1 0 0 0 0 0]))==0 && morphism(3)==0 && morphism(4)==1
    retraction1=fAtoA';
elseif sum(abs(object-[1 0 0 0 0 0]))==0 && morphism(3)==1 && morphism(4)==0
    retraction1=AftoA;
elseif sum(abs(object-[1 0 0 0 0 0]))==0 && morphism(3)==1 && morphism(4)==1
    retraction1=fAAf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% class 111-1 [0 / 0 0 / 0]
elseif sum(abs(object-[0 1 0 0 1 0]))==0
    retraction1=[1 0;0 0];
    retraction2=[0 0;0 1];

elseif sum(abs(object-[0 1 0 0 0 0]))==0 && mod(morphism(2)+morphism(4),2)==0
    retraction1=[1 0];
elseif sum(abs(object-[0 1 0 0 0 0]))==0 && mod(morphism(2)+morphism(4),2)==1
    retraction1=[0 1];

elseif sum(abs(object-[0 0 0 0 1 0]))==0 && mod(morphism(1)+morphism(3),2)==0
    retraction1=[1; 0];
elseif sum(abs(object-[0 0 0 0 1 0]))==0 && mod(morphism(1)+morphism(3),2)==1
    retraction1=[0; 1];

elseif sum(abs(object-[0 0 0 0 0 0]))==0
    retraction1=1;
end

[target_V, origin_V]=size(retraction1);

if target_V>=origin_V
    norm1=retraction1'*retraction1;
    retraction1=retraction1/sqrt(norm1(1,1));
    retraction2=retraction2/sqrt(norm1(1,1));
else
    norm1=retraction1*retraction1';
    retraction1=retraction1/sqrt(norm1(1,1));
    retraction2=retraction2/sqrt(norm1(1,1));
end


retraction1=vpa(retraction1,2);
retraction2=vpa(retraction2,2);