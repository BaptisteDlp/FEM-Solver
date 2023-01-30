#include<iostream>
#include<cmath>
#include"Penalization.hpp"

using namespace std;





/**
 *Fonction qui penalise la matrice pour des conditions de dirichlet homogene
 */
void pseudoElemA(DenseMatrix &A, vector<int> const& boundary){
   
  for(unsigned int j = 0; j<boundary.size(); j++){
        
    for(int i = 0; i<A.getnbR(); i++){
            
            
      if(boundary[j]==i){
	A(i,i) = 1;
      }else{
	A(i,boundary[j]) = 0;
	A(boundary[j],i) = 0;
      }
				
    }
  }
}






/**
*Fonction qui penalise le vecteur F du systeme lineaire
*/
void penalizationF(Vecteur &F, vector<int> const& boundary){
    
    for(unsigned int j = 0; j<boundary.size(); j++){
        F(boundary[j]) = 0;
	}
	
			
}		
			
			
	
	
	
	

/**
*Fonction qui penalise la matrice A du systeme lineaire par la methode de
* la tgv, boundary est le tableau qui contient les noeuds du bords
*/
void penalizationA(DenseMatrix &A, vector<int> const& boundary){
    
    for(unsigned int j = 0; j<boundary.size(); j++){
        A(boundary[j],boundary[j]) = pow(10,30);	
	}		
			
			
			
}
			


			
			
			
/**
*Fonction qui penalise le second membre par la methode de la tgv pour un maillage quelconque
* pour des conditions de dirichlet non homogenes, boundary est le tableau qui contient les noeuds du bords
*/			
void penalizationFNonHomo(Vecteur &F, vector<int> const& boundary,DenseMatrix const& tabNodes){
   
   for(unsigned int j = 0; j<boundary.size(); j++){
	   double x = tabNodes(boundary[j],0);
	   double y = tabNodes(boundary[j],1);
	   double h = exp(x); //Change function here
       F(boundary[j]) = h*pow(10,30); 
   }
       
}











DenseMatrix periodicConditionA(DenseMatrix const& A, vector<int> const& tabIdBord){

    DenseMatrix B(A);
	double psi(0);

    
    for(unsigned int i = 0; i<tabIdBord.size(); i++){ 
        int si = tabIdBord[i]; 
        if(si != -1){
            //si appartient a gamma0 et i appartient a gamma1 donc y_si = y_i
            for(int k = 0; k<A.getnbR(); k++){
                //on remplace la ligne et la colonne correspondant au noeud i de gamma1
                //par la contribution psi
                int l = tabIdBord[k];
                
                
                if(l != -1){
                    //l appartient a gamma0 et k appartient a gamma1
                    psi = A(si,l) + A(si,k) + A(i,l) + A(i,k);
                }else{
                    //l n appartient pas a gamma0 donc k n appartient pas a gamma1
                    psi = A(si,k) + A(i,k);
				}
                B(i,k) = psi;
                B(k,i) = psi;
				
			}
		}
	}

    for(int i = 0; i<A.getnbR(); i++){
        //on supprime la ligne et la colonne correspondant aux noeuds si de gamma0
        int si = tabIdBord[i];
        if(si != -1){
            for(int k = 0; k<A.getnbR(); k++){
                if(si != k){
                    B(si,k) = 0;
                    B(k,si) = 0;
                }else{
                    B(si,si) = 1;
				}
			}
		}
                    
	}     
	
    return B;


}





void periodicConditionA2(DenseMatrix &A, vector<int> const& tabIdBord){

	double D(0);
	double psi(0);

    for(unsigned int i = 0; i<tabIdBord.size(); i++){ 
        int si = tabIdBord[i]; 
        if(si != -1){
            //si appartient a gamma0 et i appartient a gamma1 donc y_si = y_i
            D = A(i,i) + A(si,si) +2*A(si,i);
            for(int k = 0; k<A.getnbR(); k++){
                //on remplace la ligne et la colonne correspondant au noeud i de gamma1
                //par la contribution psi
                if(i != k){
                    psi = A(i,k) + A(si,k);
                    A(i,k) = psi;
                    A(k,i) = psi;
                }else{
                    A(i,i) = D;
				}
			}
            for(int k = 0; k<A.getnbR(); k++){
                //on supprime la ligne et la colonne correspondant au noeud si de gamma0
                if(si != k){
                    A(si,k) = 0;
                    A(k,si) = 0;
                }else{
                     A(si,si) = 1;
				}
					 
			}
                     
		}
	}
}
    






void periodicConditionF(Vecteur &F, vector<int> const& tabIdBord){

     for(unsigned int i = 0; i<tabIdBord.size(); i++){ 
        int si = tabIdBord[i]; 
        if(si != -1){
            //si appartient a gamma0 et i appartient a gamma1
            //on affecte les valeurs des noeuds de gamma1 a ceux de gamma 0
            F(i) = F(i) + F(si);
            F(si) = 0; 
		}
            
	 }
     

}






void periodicConditionSol(Vecteur &Sol, vector<int> const& tabIdBord){



	for(unsigned int i = 0; i<tabIdBord.size(); i++){ 
        int si = tabIdBord[i]; 
        if(si != -1){
            //si appartient a gamma0 et i appartient a gamma1
            //on affecte les valeurs des noeuds de gamma1 a ceux de gamma 0
            Sol(si) = Sol(i);
		}
            
	}


}		
			
			
			
			
			
			
			
			
			
