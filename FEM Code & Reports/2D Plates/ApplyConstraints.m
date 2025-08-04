%Coded on Octave
%This is a generalized script that makes the necessary adjustments and equations to apply constraints. It can be used for any number of dimensions and apply any number and type of linear constraints.
function [Equations Knew Fnew]=ApplyConstraints(K,F,BCs)
  %Decoding data and creating a linear system that contains the constraint equations
  FreeDegs=rows(BCs{1})-1;
  %Building linear system
  for i=1:length(BCs)
    BC=BCs{i};
    [Waste ColPlaces1 ColPlaces2]=intersect(BC(1,1:end-1),1:columns(K)/FreeDegs);
    Const(i,1)=sum(BC(:,end));
    for j=1:FreeDegs
      LinSys(i,FreeDegs*ColPlaces2-(FreeDegs-j))=BC(j+1,ColPlaces1);
    endfor   
  endfor 
  

  %Keeping only the constrained degrees of freedom and adding the constants in the linear system
  ConstDegs=find(sum(abs(LinSys)));
  LinSys=LinSys(:,ConstDegs);
  LinSys(:,end+1)=Const;

  %Computing the reduced row echelon form, to make BCs directly applicable on the stiffness matrix
  [LinSys SlaveCols]=rref(LinSys);

  %In case of zero rows the system needs further examination and changes
  zero_rows=find(sum(transpose(abs(LinSys(:,1:end-1))))==0);
  
  if length(zero_rows)>0
    %If the system has redundant rows, those get ignored, else if the system is inconsistent, an error appears and the function seizes
    if zeros(zero_rows,1)==LinSys(zero_rows,end);
      LinSys=LinSys(setdiff([1:rows(LinSys)],zero_rows),:);
    else
      Knew=[];Fnew=[];Equations={};
      errordlg('The input constraint equations are mathematicaly inconsistent. The program has stopped.','Inconsistent Constraints')
      return
    endif
  endif
  

  %Applying the BCs with the Master-Slave method
  AllDegs=[1:rows(K)];
  SlaveDegs=ConstDegs(SlaveCols);
  MasterDegs=setdiff(ConstDegs,SlaveDegs);
  IndepDegs=setdiff(AllDegs,ConstDegs);
  %Particion of the linear system into the slave and master matrices
  Am=transpose(LinSys(:,setdiff([1:(columns(LinSys)-1)],SlaveCols)));
  G=zeros(length(AllDegs),1);
  G(SlaveDegs,1)=LinSys(:,end);
  T(union(MasterDegs,IndepDegs),union(MasterDegs,IndepDegs))=eye(length(AllDegs)-length(SlaveDegs));
  T(SlaveDegs,union(MasterDegs,IndepDegs))=zeros(length(SlaveDegs),length(AllDegs)-length(SlaveDegs));
  T(SlaveDegs,MasterDegs)=-Am;
  T=T(:,union(MasterDegs,IndepDegs));

  %Outputs
  Equations={T,G};

  Knew=transpose(T)*K*T;

  Fnew=transpose(T)*(F-K*G);

endfunction
