%Solver
%Antonios-Ioakeim Kyriakopoulos 6337
%Coded on Octave
%In the Solver, stiffness matrices are constructed. The matrices are used to solve the displacement equations.
%Then the external forces are calculated and the stresses of the bar structure. The stress calculation of the beam structure will also be added, soon? (Perhaps add Von Misses stresses too)

load PPOutput6337.txt
BarNodes;BeamNodes;Beams;Bars;BarPCs;BeamPCs;BarForces;BeamForces;

%Bar Structure Solver
%Preparing Equation Matrices
BarStiffnessMatrix=zeros(3*rows(BarNodes));
BarForceMatrix=zeros(3*rows(BarNodes),1);

%Adding the elements
for i=1:rows(Bars)
  Node1=Bars(i,1);Node2=Bars(i,2);
  E=Bars(i,3);A=Bars(i,4);
  L=norm(BarNodes(Node1,:)-BarNodes(Node2,:));
  lb=(BarNodes(Node1,1)-BarNodes(Node2,1))/L;
  mb=(BarNodes(Node1,2)-BarNodes(Node2,2))/L;
  nb=(BarNodes(Node1,3)-BarNodes(Node2,3))/L;
  B=[lb mb nb 0 0 0;0 0 0 lb mb nb];
  D=A*E/L*[1 -1;-1 1];
  K=transpose(B)*D*B;
  BarStiffnessMatrix((3*Node1-2):(3*Node1),(3*Node1-2):(3*Node1))=BarStiffnessMatrix((3*Node1-2):(3*Node1),(3*Node1-2):(3*Node1))+K(1:3,1:3);
  BarStiffnessMatrix((3*Node2-2):(3*Node2),(3*Node2-2):(3*Node2))=BarStiffnessMatrix((3*Node2-2):(3*Node2),(3*Node2-2):(3*Node2))+K(4:6,4:6);
  BarStiffnessMatrix((3*Node1-2):(3*Node1),(3*Node2-2):(3*Node2))=BarStiffnessMatrix((3*Node1-2):(3*Node1),(3*Node2-2):(3*Node2))+K(1:3,4:6);
  BarStiffnessMatrix((3*Node2-2):(3*Node2),(3*Node1-2):(3*Node1))=BarStiffnessMatrix((3*Node2-2):(3*Node2),(3*Node1-2):(3*Node1))+K(4:6,1:3); 
  BBar{i}=B;
  DBar{i}=D;
endfor

%External Forces
for i=1:rows(BarForces)
  Node1=BarForces(i,1);
  BarForceMatrix(3*Node1-2:3*Node1)=BarForceMatrix(3*Node1-2:3*Node1)+transpose(BarForces(i,2:4));
endfor

%Adding Constraints
[Equations BarStiffnessMatrix_adjusted BarForceMatrix_adjusted]=ApplyConstraints(BarStiffnessMatrix,BarForceMatrix,BarPCs);
%ApplyConstraints is a generalized function that makes the necessary adjustments and equations to apply constraints. 
%It can be used for any number of dimensions and apply any number and type of linear constraints. It is coded in a different file.

%Solving Equations
BarDisplacementMatrix_adjusted=BarStiffnessMatrix_adjusted\BarForceMatrix_adjusted;

%Computing the full displacement and force matrices
TransfMatrix=Equations{1};ConstantsMatrix=Equations{2};
BarDisplacementMatrix=TransfMatrix*BarDisplacementMatrix_adjusted+ConstantsMatrix;
BarForceMatrix=BarStiffnessMatrix*BarDisplacementMatrix;

%Rearanging Displacements and External Forces
for i=1:3
  BarDisplacements(:,i)=BarDisplacementMatrix(i:3:end);
endfor
BarForces=[];
for i=1:3
  BarForces(:,i+1)=BarForceMatrix(i:3:end);
endfor
BarForces(:,1)=transpose(1:rows(BarForces));
%BarDisplacements%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bar Stresses
BarStress=zeros(rows(Bars),1);
for i=1:rows(Bars)
  Node1=Bars(i,1);Node2=Bars(i,2);A=Bars(i,4);
  BarStress(i,1)=(DBar{i}*BBar{i}*transpose([BarDisplacements(Node1,1:3) BarDisplacements(Node2,1:3)]))(1)/A;
endfor
%Beam Structure Solver
%Preparing Equation Matrices
BeamStiffnessMatrix=zeros(6*rows(BeamNodes));
BeamForceMatrix=zeros(6*rows(BeamNodes),1);

%Adding the elements
for i=1:rows(Beams);
  Node1=Beams(i,1);Node2=Beams(i,2);
  L=norm(BeamNodes(Node2,:)-BeamNodes(Node1,:));
  RotVector=Beams(i,3:5);E=Beams(i,6);G=Beams(i,7);
  Iz=Beams(i,8);Iy=Beams(i,9);A=Beams(i,10);
  J=Iz+Iy;
  %Preparing element sub matrix
  D=[ A*E/L 0 0 0 0 0 -A*E/L 0 0 0 0 0
     0 12*E*Iz/L^3 0 0 0 6*E*Iz/L^2 0 -12*E*Iz/L^3 0 0 0 6*E*Iz/L^2
     0 0 12*E*Iy/L^3 0 -6*E*Iy/L^2 0 0 0 -12*E*Iy/L^3 0 -6*E*Iy/L^2 0
     0 0 0 G*J/L 0 0 0 0 0 -G*J/L 0 0
     0 0 -6*E*Iy/L^2 0 4*E*Iy/L 0 0 0 6*E*Iy/L^2 0 2*E*Iy/L 0
     0 6*E*Iz/L^2 0 0 0 4*E*Iz/L 0 -6*E*Iz/L^2 0 0 0 2*E*Iz/L
     -A*E/L 0 0 0 0 0 A*E/L 0 0 0 0 0
     0 -12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2 0 12*E*Iz/L^3 0 0 0 -6*E*Iz/L^2
     0 0 -12*E*Iy/L^3 0 6*E*Iy/L^2 0 0 0 12*E*Iy/L^3 0 6*E*Iy/L^2 0
     0 0 0 -G*J/L 0 0 0 0 0 G*J/L 0 0
     0 0 -6*E*Iy/L^2 0 2*E*Iy/L 0 0 0 6*E*Iy/L^2 0 4*E*Iy/L 0
     0 6*E*Iz/L^2 0 0 0 2*E*Iz/L 0 -6*E*Iz/L^2 0 0 0 4*E*Iz/L];
     
  %Preparing recoordination matrix
  x1=BeamNodes(Node1,1);x2=BeamNodes(Node2,1);x3=RotVector(1)+BeamNodes(Node1,1);
  y1=BeamNodes(Node1,2);y2=BeamNodes(Node2,2);y3=RotVector(2)+BeamNodes(Node1,2);
  z1=BeamNodes(Node1,3);z2=BeamNodes(Node2,3);z3=RotVector(3)+BeamNodes(Node1,3);
  x21=x2-x1;x31=x3-x1;y21=y2-y1;y31=y3-y1;z21=z2-z1;z31=z3-z1;
  A123=sqrt((y21*z31-y31*z21)^2+(z21*x31-z31*x21)^2+(x21*y31-x31*y21)^2);
  lx=x21/L;mx=y21/L;nx=z21/L;
  lz=(y21*z31-y31*z21)/A123;mz=(z21*x31-z31*x21)/A123;nz=(x21*y31-x31*y21)/A123;  
  ly=mz*nx-nz*mx;my=nz*lx-lz*nx;ny=lz*mx-mz*lx;
  T=[lx mx nx;ly my ny;lz mz nz];
  B=blkdiag(T,T,T,T);
  BeamRecoordSubMatrix=transpose(B)*D*B;
  %Adding element submatrix in the stiffness matrix
  BeamStiffnessMatrix((6*Node1-5):(6*Node1),(6*Node1-5):(6*Node1))=BeamStiffnessMatrix((6*Node1-5):(6*Node1),(6*Node1-5):(6*Node1))+BeamRecoordSubMatrix(1:6,1:6);
  BeamStiffnessMatrix((6*Node2-5):(6*Node2),(6*Node2-5):(6*Node2))=BeamStiffnessMatrix((6*Node2-5):(6*Node2),(6*Node2-5):(6*Node2))+BeamRecoordSubMatrix(7:12,7:12);
  BeamStiffnessMatrix((6*Node1-5):(6*Node1),(6*Node2-5):(6*Node2))=BeamStiffnessMatrix((6*Node1-5):(6*Node1),(6*Node2-5):(6*Node2))+BeamRecoordSubMatrix(1:6,7:12);
  BeamStiffnessMatrix((6*Node2-5):(6*Node2),(6*Node1-5):(6*Node1))=BeamStiffnessMatrix((6*Node2-5):(6*Node2),(6*Node1-5):(6*Node1))+BeamRecoordSubMatrix(7:12,1:6);  
  BBeam{i}=B;
  DBeam{i}=D;
endfor

%External Forces
for i=1:rows(BeamForces)
  Node1=BeamForces(i,1);
  BeamForceMatrix(6*Node1-5:6*Node1)=BeamForceMatrix(6*Node1-5:6*Node1)+transpose(BeamForces(i,2:7));
endfor

%Adding Constraints
[Equations BeamStiffnessMatrix_adjusted BeamForceMatrix_adjusted]=ApplyConstraints(BeamStiffnessMatrix,BeamForceMatrix,BeamPCs);
%ApplyConstraints is a generalized script that makes the necessary adjustments and equations to apply constraints. 
%It can be used for any number of dimensions and apply any number and type of linear constraints. It is coded in a different file.

%Solving Equations
BeamDisplacementMatrix_adjusted=BeamStiffnessMatrix_adjusted\BeamForceMatrix_adjusted;

%Computing the full displacement and force matrices
TransfMatrix=Equations{1};ConstantsMatrix=Equations{2};
BeamDisplacementMatrix=TransfMatrix*BeamDisplacementMatrix_adjusted+ConstantsMatrix;
BeamForceMatrix=BeamStiffnessMatrix*BeamDisplacementMatrix;

%Rearanging Displacements and External Forces
for i=1:6
  BeamDisplacements(:,i)=BeamDisplacementMatrix(i:6:end);
endfor
BeamForces=[];
for i=1:6
  BeamForces(:,i+1)=BeamForceMatrix(i:6:end);
endfor
BeamForces(:,1)=transpose(1:rows(BeamForces));
%BeamDisplacements%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Beam Stressses (maybe some day)
for i=1:rows(Beams)
Node1=Beams(i,1);Node2=Beams(i,2);
BeamElementForces(i,:)=(DBeam{i}*BBeam{i}*transpose([BeamDisplacements(Node1,1:6) BeamDisplacements(Node2,1:6)]));
endfor


%Outputs
save SolverOutput6337.txt BarNodes BeamNodes Beams Bars BarPCs BeamPCs BarForces BeamForces BarDisplacements BeamDisplacements BarStress

clear all 


