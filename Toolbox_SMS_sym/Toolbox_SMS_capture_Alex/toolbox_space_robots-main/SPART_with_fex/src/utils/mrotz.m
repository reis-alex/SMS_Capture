% mrotz		M = mrotz(a)		Matrice rotation d'angle a autour de Oz

function [m] = mrotz(a)

s = sin(a) ; c = cos(a);

m = [c -s 0 ; s c 0 ; 0 0 1] ;
     
  
