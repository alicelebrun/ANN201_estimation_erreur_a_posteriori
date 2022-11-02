%% =====================================================
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de la chaleur suivante stationnaire, avec conditions de
% Dirichlet non homogene
%
% | \alpha T - div(\lambda \grad T) = S,   dans \Omega=\Omega_1 U \Omega_2
% |        T = T_\Gamma,   sur le bord
%
% ou S est la source de chaleur, T_\Gamma la température exterieure
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

clear all;close all;
nom_maillage = 'geomChaleur.msh';
validation   = 'oui';
raffinement  = 'oui';
reference    = 'non';

alpha   = 1.0;
T_Gamma = 0.0;

%% lecture du maillage et affichage
%% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);
Nbpt

%% ----------------------
%% calcul des matrices EF
%% ----------------------

%% declarations
%% ------------
iKK = zeros(9*Nbtri, 1);
jKK = zeros(9*Nbtri, 1);
vKK = zeros(9*Nbtri, 1);

iMM = zeros(9*Nbtri, 1);
jMM = zeros(9*Nbtri, 1);
vMM = zeros(9*Nbtri, 1);

% boucle sur les triangles
% ------------------------
fprintf("Assemblage de la matrice, nb triangles=%d\n", Nbtri)
tic;
index = 1;
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  II=Numtri(l,:);
  S1=Coorneu(II(1),:);
  S2=Coorneu(II(2),:);
  S3=Coorneu(II(3),:);
  % calcul des matrices elementaires du triangle l 
  Kel=matK_elem(S1, S2, S3, Reftri(l));
  Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  for i=1:3
    I=II(i);
    for j=1:3
      J=II(j);
      iKK(index)=I;
      jKK(index)=J;
      vKK(index)=Kel(i,j);
      iMM(index)=I;
      jMM(index)=J;
      vMM(index)=Mel(i,j);
      index=index+1;
    end % for j
  end % for i

end % for l
KK = sparse(iKK,jKK,vKK); % matrice de rigidite
% Vérification KK: on doit trouver 0
sz=size(KK);
fprintf("Norme KK*1=%f\n", norm(KK*ones(sz(1),1)))
% return
MM = sparse(iMM,jMM,vMM); % matrice de masse
% Vérification MM: on doit retrouver 4, la surface du domaine
sz=size(MM);
fprintf("Valeur 1*MM*1=%f\n", ones(1,sz(1))*MM*ones(sz(1),1))
% return
AA = alpha*MM + KK;
toc
%% Calcul du second membre L
%% -------------------------
fprintf("Assemblage du second membre, nb noeuds=%d\n", Nbpt)
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM * FF;

%% Pseudo-élimination
%% ------------------
fprintf("Pseudo-élimination\n")
[tilde_AA, tilde_LL] = elimine(AA, LL, Refneu);
%% inversion
%% ----------
fprintf("Résolution")
spparms ('spumoni',3)
UU = tilde_AA\tilde_LL;
TT = UU+T_Gamma;


%% visualisation
%% -------------
fprintf("Affichage\n")
figure(1); affiche(TT, Numtri, Coorneu, sprintf('Uh - %s', nom_maillage));

%% Calcul de grad(u_h) par triangle
%% --------------------------------
grad_uh = zeros(Nbtri, 2);
for l=1:Nbtri

  %% Indices globaux des sommets du triangle
  II = Numtri(l,:);

  %% Coordonnees des sommets du triangle
  S1=Coorneu(II(1),:);
  S2=Coorneu(II(2),:);
  S3=Coorneu(II(3),:);
  gradS = grad_elem(S1, S2, S3);
  grad_uh(l,:) =UU(II)'*gradS;
end


%% Calcul de sigma_h comme une moyenne de -grad(u_h) par patch
%% -----------------------------------------------------------
Sigma = zeros(Nbpt,2);
Sigma(:,1) = moyenne_par_patch(-grad_uh(:,1), Numtri, Coorneu);
Sigma(:,2) = moyenne_par_patch(-grad_uh(:,2), Numtri, Coorneu);

figure(2); affiche(Sigma(:,1), Numtri, Coorneu, 'sigma_x');
figure(3); affiche(Sigma(:,2), Numtri, Coorneu, 'sigma_y');


%% Estimateur d'erreur
%% -------------------
Eta = zeros(Nbtri, 1);
Eta_r2 = 0;  %% Composante en résidu
Eta_f2 = 0;  %% Composante en flux
for l = 1:Nbtri
  %% Indices des sommets du triangle
  II = Numtri(l,:);

  %% Coordonnees des sommets du triangle
  S1=Coorneu(II(1),:);
  S2=Coorneu(II(2),:);
  S3=Coorneu(II(3),:);

  X = Coorneu(II,1);
  Y = Coorneu(II,2);

  %% Matrices élémentaires
  M_el = matM_elem(S1, S2, S3);
  g_el = grad_elem(S1, S2, S3);

  %% Composante en flux :
  %% E_f  = grad(u_h) + sigma_h
  %% E_f2 = || Ef ||^2_{L2}
  %% Le champ grad_uh est défini par triangle,
  %% alors que le champ Sigma est défini aux noeuds
  %% Pour un triangle Tl on compare donc une fonction
  %% constante et une fonction linéaire
  E_fx = grad_uh(l, 1)+Sigma(II,1);
  E_fy = grad_uh(l, 2)+Sigma(II,2);

  E_f2 = E_fx'*M_el*E_fx+E_fy'*M_el*E_fy;

  %% Composante en résidu :
  %% E_r  = f - div(sigma_h) - u_h
  %% E_r2 = || E_r ||^2_{L2}
  %% Sigma est un champ vectoriel P1
  %% Sigma(x,y)=[Sigma_x(S1) Sigma_y(S1)]w1(x,y)+
  %%            [Sigma_x(S2) Sigma_y(S2)]w2(x,y)+
  %%            [Sigma_x(S3) Sigma_y(S3)]w3(x,y)
  %% 
  gradSigmaX = g_el(:,1).*Sigma(II,1);
  gradSigmaY = g_el(:,2).*Sigma(II,2);
  divSigma   = gradSigmaX + gradSigmaY;

  E_r  = FF(II,1)-divSigma-UU(II,1);
  E_r2 = E_r'*M_el*E_r;

  %% Estimateur global
  Eta_f2 = Eta_f2 + E_f2;
  Eta_r2 = Eta_r2 + E_r2;
  Eta(l) = sqrt(E_f2) + 2*sqrt(2)/pi * sqrt(E_r2);
end

%% Estimateur global
Eta_tot = sqrt(Eta_f2) + 2*sqrt(2)/pi * sqrt(Eta_r2)

%% Lissage de l'estimateur par triangle, pour visualisation et post-traitement
Eta_av = moyenne_par_patch(Eta, Numtri, Coorneu);
figure(4); affiche(Eta_av, Numtri, Coorneu, 'Eta');


%% validation
%% ----------
if strcmp(validation,'oui')
  X = Coorneu(:,1);
  Y = Coorneu(:,2);

  %% CAS I : REFERENCE ANALYTIQUE
  UU_ref = sin(pi*X) .* sin(2*pi*Y);
  Sigma_ref = zeros(Nbpt, 2);
  Sigma_ref(:,1) = -pi*cos(pi*X).*sin(2*pi*Y); %% A COMPLETER
  Sigma_ref(:,2) = -2*pi*sin(pi*X).*cos(2*pi*Y); %% A COMPLETER

  %% CAS II : REFERENCE NON ANALYTIQUE
  %% Récupération de la référence
  %% et interpolation sur le maillage courant
  if strcmp(reference,'oui')
    [UU_ref, Sigma_ref] = lecture_reference(X, Y);
  end
  %% Erreur sur u
  EE = (UU-UU_ref);
  figure(5); affiche(abs(EE), Numtri, Coorneu, sprintf('Erreur - %s', nom_maillage));

  %% Calcul de l erreur L2 sur u_h
  err_u_l2 = EE'*MM*EE
  %% Calcul de l erreur en semi-norme H1 sur u_h
  err_u_semih1 = EE'*KK*EE


  %% Erreur sur sigma
  EEx = Sigma(:,1) - Sigma_ref(:,1);
  EEy = Sigma(:,2) - Sigma_ref(:,2);

  %% Calcul de l'erreur L2 sur sigma_h
  err_sigma_l2 = sqrt(EEx'*MM*EEx+EEy'*MM*EEy)
end

%% Raffinement du maillage
if strcmp(raffinement, 'oui')
  h = 1 ./ Eta_av;     % Calcul du nouveau pas de maillage
  h = h * 0.01/min(h); % - mise à l'échelle pour que : min(h) == 0.01
  h = min(h, 0.4);     % - écrêtement pour que :       max(h) <= 0.4
  figure(6); affiche(h, Numtri, Coorneu, 'h')

  %% Génération des fichiers pour gmsh
  %%
  %% Commande pour générer un maillage raffiné :
  %%   gmsh -bgm size.pos geomCarre.geo
  write_field('size.pos', 'Mesh size', h, Numtri, Coorneu);
end

    
%% Enregistrement des résultats de référence (si maillage suffisamment raffiné)
if strcmp(reference, 'oui')
  dx = 1e-3;
  X = Coorneu(:,1);
  Y = Coorneu(:,2);
  Xrange = min(X):dx:max(X);
  Yrange = min(Y):dx:max(Y);
  [XI,YI] = meshgrid(Xrange, Yrange);
  Uref      = griddata(X,Y,UU,        XI,YI);
  SigmaXref = griddata(X,Y,Sigma(:,1),XI,YI);
  SigmaYref = griddata(X,Y,Sigma(:,2),XI,YI);
  save 'ref.mat' Xrange Yrange Uref SigmaXref SigmaYref
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
