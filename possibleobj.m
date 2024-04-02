function [object5] = possibleobj(morphism5)
%POSSIBLEOBJ 此处显示有关此函数的摘要
%   此处显示详细说明
ob01=1;  %initial state, can be chosen,32 possible
ob12=1;
ob23=1;
ob34=1;
ob45=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
if mor012==2
    ob02=mod(3-ob01-ob12,2);
else
    ob02=mod(2-ob01-ob12,2);
end

if mor123==2
    ob13=mod(3-ob12-ob23,2);
else
    ob13=mod(2-ob12-ob23,2);
end

if mor234==2
    ob24=mod(3-ob23-ob34,2);
else
    ob24=mod(2-ob23-ob34,2);
end

if mor345==2
    ob35=mod(3-ob34-ob45,2);
else
    ob35=mod(2-ob34-ob45,2);
end

if mor023==2
    ob03=mod(3-ob02-ob23,2);
else
    ob03=mod(2-ob02-ob23,2);
end

if mor134==2
    ob14=mod(3-ob13-ob34,2);
else
    ob14=mod(2-ob13-ob34,2);
end

if mor245==2
    ob25=mod(3-ob24-ob45,2);
else
    ob25=mod(2-ob24-ob45,2);
end

if mor034==2
    ob04=mod(3-ob03-ob34,2);
else
    ob04=mod(2-ob03-ob34,2);
end

if mor145==2
    ob15=mod(3-ob14-ob45,2);
else
    ob15=mod(2-ob14-ob45,2);
end

if mor045==2
    ob05=mod(3-ob04-ob45,2);
else
    ob05=mod(2-ob04-ob45,2);
end

object5=[ob01 ob02 ob03 ob04 ob05 ob12 ob13 ob14 ob15 ob23 ob24 ob25 ob34 ob35 ob45];


end

