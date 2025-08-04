%OCTAVE
load PPOutput6337.txt
Elements;Nodes;Forces;PCs;

StiffnessMatrix=zeros(2*rows(Nodes));
ForceMatrix=zeros(2*rows(Nodes),1);

for i=1:rows(Elements)

  Node1=Elements(i,1);Node2=Elements(i,2);Node3=Elements(i,3);
  E=Elements(i,5);v=Elements(i,6);t=Elements(i,4);
  x1=Nodes(Node1,1);x2=Nodes(Node2,1);x3=Nodes(Node3,1);
  y1=Nodes(Node1,2);y2=Nodes(Node2,2);y3=Nodes(Node3,2);
  y23=y2-y3;y31=y3-y1;y12=y1-y2;y13=y1-y3;
  x32=x3-x2;x13=x1-x3;x21=x2-x1;x23=x2-x3;
  A=polyarea([x1;x2;x3],[y1;y2;y3]);
  B=(2*A)^-1*[y23 0 y31 0 y12 0;0 x32 0 x13 0 x21;x32 y23 x13 y31 x21 y12];
  %D=E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
  D=E/((1+v)*(1-2*v))*[1-v v 0;v 1-v 0;0 0 (1-2*v)/2];
  K=transpose(B)*D*B*t*A;
  StiffnessMatrix(2*Node1-1:2*Node1,2*Node1-1:2*Node1)=StiffnessMatrix(2*Node1-1:2*Node1,2*Node1-1:2*Node1)+K(1:2,1:2);
  StiffnessMatrix(2*Node1-1:2*Node1,2*Node2-1:2*Node2)=StiffnessMatrix(2*Node1-1:2*Node1,2*Node2-1:2*Node2)+K(1:2,3:4);
  StiffnessMatrix(2*Node1-1:2*Node1,2*Node3-1:2*Node3)=StiffnessMatrix(2*Node1-1:2*Node1,2*Node3-1:2*Node3)+K(1:2,5:6);
  StiffnessMatrix(2*Node2-1:2*Node2,2*Node1-1:2*Node1)=StiffnessMatrix(2*Node2-1:2*Node2,2*Node1-1:2*Node1)+K(3:4,1:2);
  StiffnessMatrix(2*Node2-1:2*Node2,2*Node2-1:2*Node2)=StiffnessMatrix(2*Node2-1:2*Node2,2*Node2-1:2*Node2)+K(3:4,3:4);
  StiffnessMatrix(2*Node2-1:2*Node2,2*Node3-1:2*Node3)=StiffnessMatrix(2*Node2-1:2*Node2,2*Node3-1:2*Node3)+K(3:4,5:6);
  StiffnessMatrix(2*Node3-1:2*Node3,2*Node1-1:2*Node1)=StiffnessMatrix(2*Node3-1:2*Node3,2*Node1-1:2*Node1)+K(5:6,1:2);
  StiffnessMatrix(2*Node3-1:2*Node3,2*Node2-1:2*Node2)=StiffnessMatrix(2*Node3-1:2*Node3,2*Node2-1:2*Node2)+K(5:6,3:4);
  StiffnessMatrix(2*Node3-1:2*Node3,2*Node3-1:2*Node3)=StiffnessMatrix(2*Node3-1:2*Node3,2*Node3-1:2*Node3)+K(5:6,5:6);
  Dcell{i}=D;
  Bcell{i}=B;
endfor
%External Forces
for i=1:rows(Forces)
  Node1=Forces(i,1);
  ForceMatrix(2*Node1-1:2*Node1)=ForceMatrix(2*Node1-1:2*Node1)+transpose(Forces(i,2:3));
endfor
%Constraints
[Equations StiffnessMatrix_adjusted ForceMatrix_adjusted]=ApplyConstraints(StiffnessMatrix,ForceMatrix,PCs);
%Solver

DisplacementMatrix_adjusted=StiffnessMatrix_adjusted\ForceMatrix_adjusted;
%Compute the full displacement and force matrices

TransfMatrix=Equations{1};ConstantsMatrix=Equations{2};
DisplacementMatrix=TransfMatrix*DisplacementMatrix_adjusted+ConstantsMatrix;

ForceMatrix=StiffnessMatrix*DisplacementMatrix;


%Transforming Displacements
for i=1:2
  Displacements(:,i)=DisplacementMatrix(i:2:end);
endfor
%Calculating Stresses and Deformations

for i=1:rows(Elements)

  Node1=Elements(i,1);
  Node2=Elements(i,2);
  Node3=Elements(i,3);
  ElementStress(1:3)=transpose(Dcell{i}*Bcell{i}*transpose([Displacements(Node1,:) Displacements(Node2,:) Displacements(Node3,:)]));
  Deformation(1:3)=transpose(Bcell{i}*transpose([Displacements(Node1,:) Displacements(Node2,:) Displacements(Node3,:)]));
  ex(i,1)=Deformation(1);ey(i,1)=Deformation(2);gxy(i,1)=Deformation(3);sx(i,1)=ElementStress(1);sy(i,1)=ElementStress(2);txy(i,1)=ElementStress(3);
  Svm(i,1)=sqrt((sx(i)^2)+(sy(i)^2)-(sx(i))*(sy(i))+3*(txy(i)^2));
  s1(i,1)=((sx(i)+sy(i))/2)+sqrt((((sx(i)-sy(i))/2)^2)+txy(i)^2);
  s2(i,1)=((sx(i)+sy(i))/2)-sqrt((((sx(i)-sy(i))/2)^2)+txy(i)^2);
  thita(i,:)=atand(2*txy(i)/(sx(i)-sy(i)))/2;
endfor

Deformations={ex ey gxy};
Stresses={Svm sx sy txy};
PrincStresses={s1 s2};
save SolverOutput6337.txt Elements Nodes Forces PCs Displacements Deformations Stresses PrincStresses thita

clear all











