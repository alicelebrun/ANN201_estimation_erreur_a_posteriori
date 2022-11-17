function [h_min, h_max]=calcul_h(Nbtri, Coorneu, Numtri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul_h.m:
% routine de calcul du rayon maximum des cercles circonscrits aux triangles
% d'un maillage
%
% SYNOPSIS [h_min, h_max]=calcul_h(Nbtri, Numtri)
%
% INPUT * Nbtri : nbre de triangles (entier)
%       * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%       * Numtri : liste de triangles
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%
% OUTPUT - h_min : plus petit rayon de cercle circonscrit des triangles
%        - h_max : plus grand rayon de cercle circonscrit des triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
h_min=inf;
h_max=0.0;
for l=1:Nbtri
  II=Numtri(l,:);
  S1=Coorneu(II(1),:);
  S2=Coorneu(II(2),:);
  S3=Coorneu(II(3),:);
  a=norm(S1-S2);
  b=norm(S2-S3);
  c=norm(S3-S1);
  p=0.5*(a+b+c);
  r=a*b*c/(4*sqrt(p*(p-a)*(p-b)*(p-c)));
  h_min=min(h_min, r);
  h_max=max(h_max, r);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
