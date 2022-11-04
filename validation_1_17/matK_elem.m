function [Kel] = matK_elem(S1, S2, S3, Reftri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrice de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul par passage a l'element de reference
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);
% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);

% Points et poids de quadrature
M_hat(:,1) = [1/6; 1/6];
M_hat(:,2) = [2/3; 1/6];
M_hat(:,3) = [1/6, 2/3];
poids = 1/6;

% Matrice de transport de l'élément de référence vers
% l'élément du maillage
B_l = [S2 - S1;S3 - S1];
% Gradient des fonctions de base
gradW = grad_elem(S1, S2, S3)';
gradWIgradWJ = poids*D*gradW'*gradW;
Kel = zeros(3,3);
for i=1:3
  for j=1:3
    dotGrad = gradWIgradWJ(i,j);
    % On transporte en une fois tous les points d'integration
    M_tri = B_l*M_hat+S1';
    if Reftri == 1
      lambdaNode = lambda1(M_tri(1,:), M_tri(2,:));
    else
      lambdaNode = lambda2(M_tri(1,:), M_tri(2,:));
    end
    for k=1:3
      % Evaluation du coefficient de diffusion.
      Kel(i,j) = Kel(i,j)+lambdaNode(k)*dotGrad;
    end % k
  end % j
end % i


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
