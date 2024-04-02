function [interchange] = interchange(object4,morphism4)
ob01=object4(1);
ob02=object4(2);
%ob03=object4(3);
ob04=object4(4);
ob12=object4(5);
%ob13=object4(6);
%ob14=object4(7);
ob23=object4(8);
ob24=object4(9);
ob34=object4(10);
mor012=morphism4(1);
%mor013=morphism4(2);
%mor014=morphism4(3);
%mor023=morphism4(4);
%mor024=morphism4(5);
%mor034=morphism4(6);
%mor123=morphism4(7);
%mor124=morphism4(8);
%mor134=morphism4(9);
mor234=morphism4(10);

%check the dimension of the bimodule (after quotient)
%for 012 V1
if ob01+ob02+ob12==3 %AA 1
    V1=2;
elseif ob01+ob02+ob12==2 && ob02==0 % V 1
    V1=2;
elseif ob01+ob02+ob12==1 && ob02==0 % _A A _1 2
    V1=2;
else
    V1=1; % _A A _A 2; _1 A _A 1;  1
end

%for 234 V2
if ob23+ob24+ob34==3 %AA 1
    V2=2;
elseif ob23+ob24+ob34==2 && ob24==0 % V 1
    V2=2;
elseif ob23+ob24+ob34==1 && ob24==0 % _A A _1 2
    V2=2;
else
    V2=1; % _A A _A 2; _1 A _A 1;  1
end

%for 024 V3
if ob02+ob04+ob24==3 %AA 1
    V3=4;
elseif ob02+ob04+ob24==0
    V3=1;
else
    V3=2;
end
%%%%%%%%%%%%%%% 22
if V1*V2==4 && mor012+mor234==0
    upint=[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1];
elseif V1*V2==4 && mor012==0 && mor234==1
    upint=[1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
elseif V1*V2==4 && mor012==1 && mor234==0
    upint=[1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1];
elseif V1*V2==4 && mor012==1 && mor234==1
    upint=[-1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
%%%%%%%%%%%%%%% 21
elseif V1==2 && V2==1 && mor012==0 && mor234==0
    upint=eye(2);
elseif V1==2 && V2==1 && mor012==0 && mor234==1
    upint=[1 0; 0 -1];
elseif V1==2 && V2==1 && mor012==1 && mor234==0
    upint=eye(2);
elseif V1==2 && V2==1 && mor012==1 && mor234==1
    upint=[-1 0; 0 1];
%%%%%%%%%%%%%%% 12
elseif V1==1 && V2==2 && mor012==0 && mor234==0
    upint=eye(2);
elseif V1==1 && V2==2 && mor012==0 && mor234==1
    upint=eye(2);
elseif V1==1 && V2==2 && mor012==1 && mor234==0
    upint=[1 0; 0 -1];
elseif V1==1 && V2==2 && mor012==1 && mor234==1
    upint=[-1 0; 0 1];
%%%%%%%%%%%%%%% 11
elseif V1==1 && V2==1 && mor012==1 && mor234==1
    upint=-1;
else
    upint=1;
end

interchange=kron(upint,eye(V3));
end

