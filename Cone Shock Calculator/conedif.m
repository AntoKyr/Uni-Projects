function conical=conedif(theta,z)
  gamma=1.4;
  val=(gamma-1)/2;
  conical(1,1)=z(2);
  conical(2,1)=((-2*val*z(1))-(val*z(2)*cot(theta))+(2*val*z(1)^3)+(val*z(1)^2*z(2)*cot(theta))+(2*val*z(1)*z(2)^2)+(val*z(2)^3*cot(theta))+(z(1)*z(2)^2))/(val*(1-z(1)^2-z(2)^2)-z(2)^2);
endfunction
