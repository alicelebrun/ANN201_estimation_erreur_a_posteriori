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

% Extraction de la diagonale de AA assocée aux noeud Dirichlet
indices_D = find(Refneu==1);
LLE(indices_D) = 0.0;
n = size(AA,1);
fprintf("Elimination de %d noeuds\n", size(indices_D,1));
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
