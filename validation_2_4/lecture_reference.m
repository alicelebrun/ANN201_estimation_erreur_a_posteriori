function [UU_ref, Sigma_ref] = lecture_reference(X, Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lecture_reference : Lecture d'un fichier contenant une solution de référence
% (décrite par une grille fine dans laquelle on interpole les valeurs aux points
% du maillage courant)
%
% SYNOPSIS [UU_ref, Sigma_ref] = lecture_reference(filename, X, Y)
%
% INPUT * filename : nom du fichier '.mat' contenant la solution de référence.
%       * X, Y : coordonnées des noeuds du maillage courant.
%
% OUTPUT * UU_ref : valeur de référence de la solution (u)
%        * Sigma_ref : valeur de référence du flux associé à la solution (sigma=-grad(u))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Chargement des données du fichier
  load 'ref.mat';
  Sigma_ref = zeros(length(X), 2);

  %% Reconstruction du maillage fin
  [XI, YI]  = meshgrid(Xrange, Yrange);

  %% Interpolation des différentes quantités intéressantes
  UU_ref         = interp2(XI, YI, Uref,      X, Y);
  Sigma_ref(:,1) = interp2(XI, YI, SigmaXref, X, Y);
  Sigma_ref(:,2) = interp2(XI, YI, SigmaYref, X, Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
