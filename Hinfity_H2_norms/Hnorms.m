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
  
  G = ss(A,B,C,D)
  [GS,GNS]=stabsep(G)
  
  
  GSP = GS.B*(GS.B)';
  GSL = lyap(GS.A,GSP);
  GSnorm = trace(GS.C*GSL*GS.C');
  
  GNSP = GNS.B*(GNS.B)';
  GNSL = lyap(-GNS.A,GNSP);
  GNSnorm = trace(GNS.C*GNSL*GNS.C');
  
  g2 = sqrt(GSnorm + GNSnorm);

endfunction
