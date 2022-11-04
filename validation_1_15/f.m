function val = f(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Validation
  %% val = ones(size(x)) ;

  %% Validation (question 1.13)
  val = (1+5*pi^2)*sin(pi*x).*sin(2*pi*y);

  %% Résolution numérique (question 1.16)
  %% val = 290 + 50*(abs(x-1)<0.3) .* (abs(y-1)<0.2);

  %% Etude de l'estimateur a posteriori (question 2.7)
  %% val = 50*(abs(x-0.4)<0.2) .* (abs(y-1.4)<0.15);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022
