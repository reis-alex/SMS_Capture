% mrotx		M = mrotx(a)		Matrice rotation d'angle a autour de Ox

function [m] = mrotx(a)

s = sin(a) ; c = cos(a);

m = [1 0 0 ; 0 c -s ; 0 s c] ;
     
  
