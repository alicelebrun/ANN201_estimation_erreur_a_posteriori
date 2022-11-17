function [grad] = grad_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grad_elem :
% calcul du gradient des fonctions de base P1 dans un triangle
%
% SYNOPSIS [G] = grad_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - gradient des fonctions de base P1 restreintes au triangle
%            G(i,j) = ∂wⱼ / ∂xᵢ,    ∀ i∈1:2, ∀ j∈1:3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end;

grad = norm / D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
% Dans le plan de référence, w1(x,y)=1-x-y, w2(x,y)=x, w3(x,y)=y
% Les gradients sont donc:
% gradRef=[-1 1 0]
%         [-1 0 1]
% La matrice Bl s'écrit
% [x2-x1 x3-x1]
% [y2-y1 y3-y1]
% dont l'inverse de la transposée est
% 1/D [y3-y1 y1-y2]
%     [x1-x3 x2-x1]
% avec D=(x2-x1)(y3-y1)-(x3-x1)(y2-y1)
% On en déduit les gradients dans le plan physique
% gradPhy=Bl^{-T}gradPhys
%        =[y2-y3 y3-y1 y1-y2]
%         [x3-x2 x1-x3 x2-x1]