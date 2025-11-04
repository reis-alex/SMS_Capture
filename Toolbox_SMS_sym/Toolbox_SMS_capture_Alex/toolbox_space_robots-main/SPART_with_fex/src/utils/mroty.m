% mroty		M = mroty(a)		Matrice rotation d'angle a autour de Oy

function [m] = mroty(a)

s = sin(a) ; c = cos(a);

m = [c 0 s ; 0 1 0 ; -s 0 c] ;
     
  
