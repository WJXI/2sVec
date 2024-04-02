function [Map135] = hexagon135(object5,morphism5)
ob01=object5(1);
ob02=object5(2);
ob03=object5(3);
ob04=object5(4);
ob05=object5(5);
ob12=object5(6);
ob13=object5(7);
ob14=object5(8);
ob15=object5(9);
ob23=object5(10);
ob24=object5(11);
ob25=object5(12);
ob34=object5(13);
ob35=object5(14);
ob45=object5(15);
%
mor012=morphism5(1);
mor013=morphism5(2);
mor014=morphism5(3);
mor015=morphism5(4);
mor023=morphism5(5);
mor024=morphism5(6);
mor025=morphism5(7);
mor034=morphism5(8);
mor035=morphism5(9);
mor045=morphism5(10);
mor123=morphism5(11);
mor124=morphism5(12);
mor125=morphism5(13);
mor134=morphism5(14);
mor135=morphism5(15);
mor145=morphism5(16);
mor234=morphism5(17);
mor235=morphism5(18);
mor245=morphism5(19);
mor345=morphism5(20);


%object01234=[ob01 ob02 ob03 ob04 ob12 ob13 ob14 ob23 ob24 ob34];
%morphism01234=[mor012 mor013 mor014 mor023 mor024 mor034 mor123 mor124 mor134 mor234];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
object01235=[ob01 ob02 ob03 ob05 ob12 ob13 ob15 ob23 ob25 ob35];
morphism01235=[mor012 mor013 mor015 mor023 mor025 mor035 mor123 mor125 mor135 mor235];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%object01245=[ob01 ob02 ob04 ob05 ob12 ob14 ob15 ob24 ob25 ob45];
%morphism01245=[mor012 mor014 mor015 mor024 mor025 mor045 mor124 mor125 mor145 mor245];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
object01345=[ob01 ob03 ob04 ob05 ob13 ob14 ob15 ob34 ob35 ob45];
morphism01345=[mor013 mor014 mor015 mor034 mor035 mor045 mor134 mor135 mor145 mor345];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%object02345=[ob02 ob03 ob04 ob05 ob23 ob24 ob25 ob34 ob35 ob45];
%morphism02345=[mor023 mor024 mor025 mor034 mor035 mor045 mor234 mor235 mor245 mor345];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
object12345=[ob12 ob13 ob14 ob15 ob23 ob24 ob25 ob34 ob35 ob45];
morphism12345=[mor123 mor124 mor125 mor134 mor135 mor145 mor234 mor235 mor245 mor345];

%[tj01234,VL01234,VR01234]=tj(object01234,morphism01234); %024 1

[tj01235,VL01235,VR01235]=tj(object01235,morphism01235); %135 3

%[tj01245,VL01245,VR01245]=tj(object01245,morphism01245); %024 2

[tj01345,VL01345,VR01345]=tj(object01345,morphism01345); %135 2

%[tj02345,VL02345,VR02345]=tj(object02345,morphism02345); %024 3

[tj12345,VL12345,VR12345]=tj(object12345,morphism12345); %135 1


if VR12345(1)==2 && VR01345(2)==2
    EX1=[1 0 0 0;0 0 1 0; 0 1 0 0;0 0 0 1];
else 
    EX1=eye(VR12345(1)*VR01345(2));
end
bubble1_135=kron(kron(eye(VR01235(1)*VR01345(1)),EX1),eye(VR12345(2)*VR12345(3)));
Map1_135=kron(eye(VR01235(1)*VR01345(1)*VR01345(2)),tj12345);
Map2_135=kron(kron(eye(VR01235(1)),tj01345),eye(VL12345(2)));
if VL01345(1)==2 && VR01235(1)==2
    EX2=[1 0 0 0;0 0 1 0; 0 1 0 0;0 0 0 1];
else 
    EX2=eye(VL01345(1)*VR01235(1));
end
bubble2_135=kron(EX2,eye(VR01235(2)*VR01235(3)));
Map3_135=kron(eye(VL01345(1)),tj01235);

Map135=Map3_135*bubble2_135*Map2_135*Map1_135*bubble1_135;
Map135=vpa(Map135,2);
Map135=round(Map135,2);
end

