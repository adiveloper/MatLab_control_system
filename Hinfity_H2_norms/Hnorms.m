function [r , g2] = Hnorms(A,B,C,D)
  rL = 0.001 ; 
  rU = 50000; 
  tol = 1.0e-10 ;
  for i = 1:500
  if (rU-rL)/rL <= tol , break , end
  r= (rU+rL)/2;
  
  R = r^2*eye( size(D'*D,1) ) - D'*D ; 
  AA = A+B*inv(R)*D'*C  ;  DD = D*inv(R)*D' ;
  %H --> Hamiltonian Matrix
  H = [ AA  B*inv(R)*B'  ;  -C'*(eye(size(DD))+DD)*C  -AA' ] ;
  eigH = eig(H)' ;
  
  if min( abs(real(eigH)) ) < 1.0e-5 , rL = r;
  else rU = r;
  end
end
  
  P = B*B';
  Lc = lyap(A,P);
  g2 = sqrt(trace(C*Lc*C'));

endfunction
