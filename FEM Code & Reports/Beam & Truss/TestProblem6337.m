%This simple structure exists for the sole reason of testing the Solver against another already established Solver 

E=210*10^9;
v=0.3;
G=E/(2*(1+v));
F=10^3;

%Bar Test Problem
BarNodes=[0 0 0;1 0 0];
Bars=[1 2 E 0.05^2];
BarPCs{1}=[1 0;1 0;0 0;0 0];
BarPCs{2}=[1 0;0 0;1 0;0 0];
BarPCs{3}=[1 0;0 0;0 0;1 0];
BarPCs{4}=[2 0;0 0;0 0;1 0];
BarPCs{5}=[2 0;0 0;1 0;0 0];
BarForces=[2 100000 0 0];
%Beam Exemplar Problem
if 1==1
BeamNodes=[0 0 0;1 0 0];
Beams=[1 2 [0 0 1] E G 0.05^4/12 0.05^4/12 0.05^2];
BeamForces=[2 0 0 0 0 1000 0];
BeamPCs{1}=[1 0;1 0;0 0;0 0;0 0;0 0;0 0];
BeamPCs{2}=[1 0;0 0;1 0;0 0;0 0;0 0;0 0];
BeamPCs{3}=[1 0;0 0;0 0;1 0;0 0;0 0;0 0];
BeamPCs{4}=[1 0;0 0;0 0;0 0;1 0;0 0;0 0];
BeamPCs{5}=[1 0;0 0;0 0;0 0;0 0;1 0;0 0];
BeamPCs{6}=[1 0;0 0;0 0;0 0;0 0;0 0;1 0];
else 
BeamNodes=[transpose(linspace(0,0.5,26)) zeros(26,2)];
Beams1=[transpose([1:15;2:16]) ones(15,1)*[[0 0 1] E G 0.08^3*0.1/12 0.1^3*0.08/12 0.1*0.08]];
Beams2=[transpose([16:25;17:26]) ones(10,1)*[[0 0 1] E G 0.08^3*0.06/12 0.06^3*0.08/12 0.08*0.06]];
Beams=[Beams1;Beams2];
%First edge
BeamPCs{1}=[1 0;1 0;0 0;0 0;0 0;0 0;0 0];
BeamPCs{2}=[1 0;0 0;1 0;0 0;0 0;0 0;0 0];
BeamPCs{3}=[1 0;0 0;0 0;1 0;0 0;0 0;0 0];
BeamPCs{4}=[1 0;0 0;0 0;0 0;1 0;0 0;0 0];
BeamPCs{5}=[1 0;0 0;0 0;0 0;0 0;1 0;0 0];
BeamPCs{6}=[1 0;0 0;0 0;0 0;0 0;0 0;1 0];
%Applying forces
Load=50000;
BeamForces1=[1 0 0 0 0 Load*0.02^2/12 0;26 0 0 0 -Load*0.02^2/12 0 0;26 0 -2000 0 0 0 0;26 0 0 -Load*0.01 0 0 0];
BeamForces2=[transpose([1:25]) zeros(25,2) ones(25,1)*-Load*0.02 zeros(25,3)];
BeamForces=[BeamForces1;BeamForces2];
endif
save PPOutput6337.txt BarNodes BeamNodes Beams Bars BarPCs BeamPCs BarForces BeamForces 

clear all