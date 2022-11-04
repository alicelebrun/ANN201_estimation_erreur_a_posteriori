function [AAE,LLE] = elimine(AA,LL,Refneu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine :
% Routine qui fait la pseudo-elimination pour les noeuds avec Refneu(i) = 1
% Conditions aux limites Dirichlet homogenes.
%
% SYNOPSIS elimine(AA,LL,Refneu)
%
% INPUT * AA, LL, Refneu: La matrice et le second membre associes au
% probleme sans elimination et references des noeuds sur la frontiere
% Dirichlet.
%
% OUTPUT - AAE, LLE: On rend la matrice et le second membre une fois la
% pseudo elimination realisee.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On copie la matrice AA et le second membre dans AAE et LLE.
AAE = AA;
LLE = LL;

% La matrice AA est partitionnée en 4 blocs
% A=[A_II | A_DI]
%   [-----------]
%   [A_ID | A_DD]
% et le second membre en 2 blocs
% L=[L_I]
%   [---]
%   [L_D]
% La solution est également décomposée en 2 blocs:
% U=[U_I]
%   [---]
%   [U_D]
% La pseudo-élimination revient à imposer la valeur donnée u_D à U_D et à
% traduire cette valeur connue dans A. Ceci se fait en remplaçant A_DD par
% sa diagonale et en annulant A_ID. On obtient ainsi le système linéaire
% suivant:
% [A_II | A_DI][U_I] [   L_I  ]
% [-----------][---]=[--------]
% [  0  [ d_DD][U_D] [d_DD.u_D]
% qui se réécrit:
% [A_II |  0  ][U_I] [L_I-A_DI.u_D]
% [-----------][---]=[------------]
% [  0  [ d_DD][U_D] [  d_DD.u_D  ]
% La matrice à gauche est Ae, le second membre est Le.
%
% Dans le cas d'une condition aux limites de Dirichlet homogène, u_D=0
% donc on a:
% Le=[L_I]
%    [---]
%    [ 0 ]

% Extraction de la diagonale de AA assocée aux noeud Dirichlet
indices_D = find(Refneu==1);
LLE(indices_D) = 0.0;
n = size(AA,1);
fprintf("Elimination de %d noeuds sur %d\n", size(indices_D,1), size(Refneu, 1));
for i=indices_D
    for j=1:n
        if i ~= j
            AAE(i, j)=0.0;
            AAE(j, i)=0.0;
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
