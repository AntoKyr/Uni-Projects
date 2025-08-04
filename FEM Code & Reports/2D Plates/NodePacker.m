function Nodes=NodePacker(k)
  for i=1:2:columns(k)
    NodePackage((i-1)*rows(k)/2+1:(i+1)*rows(k)/2,1:2)=k(:,i:i+1);
  endfor
  Nodes=NodePackage;
endfunction 