# ANN201_estimation_erreur_a_posteriori
TP2 ANN201
Ce TP porte sur la programmation Matlab d'un logiciel de résolution de l'équation de la chaleur en régime stationnaire en dimension 2 à l'aide d'éléments finis $P^1$. Il comporte deux parties:
  1. La programmation du logiciel de résolution pour l'équation $\alpha T-\mathrm{div}(\lambda\nabla T)=S$ dans $\Omega$ avec $T=T_\Gamma$ sur $\partial\Omega$, où $\alpha>0$ est une constante et $\lambda:\Omega\rightarrow \mathbb{R}$ est une fonction **régulière par morceaux** vérifiant $\exists \lambda_{min},\lambda_{max}>0$ tels que $\forall x\in\Omega$, $\lambda_{min}\leq\lambda(x)\leq\lambda_{max}$. Les données sont la source $S\in L^2
\Omega)$ et $T_\Gamma\in\mathbb{R}$.
     On peut voir la solution convergée pour un domaine rectangulaire avec une source discontinue gaussienne sur le domaine privé d'un sous-rectangle et constante sur le sous-rectangle:
     <img width="874" height="656" alt="convergence_1_17_6" src="https://github.com/user-attachments/assets/b2fbf32f-65b7-4278-a293-c19bcb0a5980" />
  2. L'implémentation d'un estimateur d'erreur à postériori permettant une adaptation du maillage afin d'atteindre une précision donnée à moindre coût de calcul.
     On voit sur les figures suivantes le maillage initial et le maillage raffiné permettant d'obtenir une erreur de $10^{-2}$:
     <img width="874" height="656" alt="maillage_2_10" src="https://github.com/user-attachments/assets/707672cf-bcd2-4921-9d49-79cf45db82bb" />
     <img width="874" height="656" alt="maillage_raffine_2_10" src="https://github.com/user-attachments/assets/10d2e12c-0516-48ed-98aa-90f9e68d960a" />

     L'intérêt du raffinement se voit sur les courbes de convergence suivantes:
     <img width="874" height="656" alt="convergence_2_11" src="https://github.com/user-attachments/assets/03aa24ee-3522-4de4-8e8c-60890b87f6ae" />

     Avec un maillage uniforme à 2000 noeuds on atteint une erreur $L^2$ de 0.2 alors qu'avec un maillage adapté comportant le même nombre de noeuds on atteint une erreur $L^2$ de 0.02, soit un gain d'un ordre de grandeur.

     

