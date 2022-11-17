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
list_suffixes = ["2";"1";"05";"02";"01";"005"];
iSmax=size(list_suffixes,2)
erreurL2List = zeros(iSmax, 1);
erreurH1List = zeros(iSmax, 1);
hMinInverseList = zeros(iSmax, 1);
hMaxInverseList = zeros(iSmax, 1);
alpha   = 1.0;
T_Gamma = 0.0;
for iS=1:iSmax
    fprintf("######################\n");
    nom_maillage = strcat('../maillages/geomChaleur0', list_suffixes(iS), '.msh');

    %% lecture du maillage et affichage
    %% ---------------------------------
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);
    [h_min, h_max]=calcul_h(Nbtri, Coorneu, Numtri);
    hMinInverseList(iS)=1.0/h_min;
    hMaxInverseList(iS)=1.0/h_max;
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
    MM = sparse(iMM,jMM,vMM); % matrice de masse
    AA = alpha*MM + KK;
    %% Calcul du second membre L
    %% -------------------------
    FF = f(Coorneu(:,1),Coorneu(:,2));
    LL = MM * FF;

    %% Pseudo-élimination
    %% ------------------
    [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu);
    %% inversion
    %% ----------
    %%spparms ('spumoni',0)
    UU = tilde_AA\tilde_LL;
    t=toc;
    fprintf("t=%.3gs\n", t)    %% validation
    clf
    affiche(UU, Numtri, Coorneu, sprintf('U_h - %s', nom_maillage));
    print(gcf,strcat('U_1_15_', num2str(iS)), '-dpdf')
    print(gcf,strcat('U_1_15_', num2str(iS)), '-dpng')
    %% ----------
    X = Coorneu(:,1);
    Y = Coorneu(:,2);

    %% CAS I : REFERENCE ANALYTIQUE
    UU_ref = sin(pi*X) .* sin(2*pi*Y);

    %% Erreur sur u
    EE = (UU-UU_ref);

    %% Calcul de l erreur L2 sur u_h
    erreurL2List(iS, 1) = sqrt(EE'*MM*EE);
    %% Calcul de l erreur en semi-norme H1 sur u_h
    erreurH1List(iS, 1) = sqrt(EE'*KK*EE);
    fprintf("h_min=%.3e h_max=%.3e", h_min, h_max)
    fprintf(" Erreur L2=%.3e Erreur semi-norme H1=%.3e\n", ...
        erreurL2List(iS, 1), erreurH1List(iS, 1))
end

format long e
clf
figure
x = log10(hMaxInverseList);
yL2 = log10(erreurL2List);
pL2 = polyfit(x, yL2, 1);
plot(x, yL2, "*-r");
hold on
plot(x, polyval(pL2, x), "--r")
yH1 = log10(erreurH1List);
pH1 = polyfit(x, yH1, 1);
plot(x, yH1, "o-b");
plot(x, polyval(pH1, x), "--b");
grid on
legend('Norme L2',strcat('Norme L2 lin m=',num2str(pL2(1))),...
'Semi-norme H1',strcat('Semi-norme H1 lin m=',num2str(pH1(1))))
xlabel("log10(1/h)")
ylabel("log10(erreur)")
print(gcf,'convergence_1_15', '-dpdf')
print(gcf,'convergence_1_15', '-dpng')
save all_1_15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
