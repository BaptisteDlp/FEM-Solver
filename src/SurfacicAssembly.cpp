#include<iostream>
#include<cmath>
#include"SurfacicAssembly.hpp"


using namespace std;


/**
*Fonction qui calcule la longueur d une arete
*/
double lengthEdge(double node1X, double node1Y, double node2X , double node2Y){
	
    return sqrt(pow(node2X-node1X,2)+pow(node2Y-node1Y,2));
}




/**
*Fonction qui retourn les valeurs de la quadrature elementaire g*phi sur
*l arete constituee de node1 et node2
*/
Vecteur boundElemLin(double node1X, double node1Y, double node2X , double node2Y){
    
    Vecteur Le(2);
  

    double x = 0.5*(node1X+node2X);
    double y = 0.5*(node1Y+node2Y);
    double g = x*y; //changer fonction ici
    
    Le(0) = 0.5*g*lengthEdge(node1X,node1Y,node2X,node2Y);
	Le(1) = 0.5*g*lengthEdge(node1X,node1Y,node2X,node2Y);

    return Le;

}





/**
*Fonction qui assemble le vecteur des termes surfaciques g*phi de la forme lineaire
*dans le cas par exemple de conditions de neumann non homogene
*/
Vecteur assemblyL(DenseMatrix const& tabNodes,vector<Vecteur> const& tabEdges){
   
   
    
    Vecteur L(tabNodes.getnbR());
   


    for(unsigned int i = 0; i<tabEdges.size(); i++){ 
	
        double node1X = tabNodes(tabEdges[i](0),0);
		double node1Y = tabNodes(tabEdges[i](0),1);
        double node2X = tabNodes(tabEdges[i](1),0);
		double node2Y = tabNodes(tabEdges[i](1),1);
        
        Vecteur Le = boundElemLin(node1X,node1Y,node2X,node2Y);
 
        L(tabEdges[i](0)) += Le(0);
		L(tabEdges[i](1)) += Le(1);

			
	}
	
	
    return L;


}





/**
*Fonction qui retourne les valeurs de l'integrale surfacique phi*phi sur l arete
*constituee de node1 et node2
*/
DenseMatrix boundElemBil(double node1X, double node1Y, double node2X , double node2Y){
    
    DenseMatrix Be(2,2);

    for(int i = 0; i<2; i++){
        for(int j = 0; j<2; j++){
            if(i == j){
                Be(i,j) = lengthEdge(node1X,node1Y,node2X,node2Y)/4;
            }else{
                Be(i,j) =  lengthEdge(node1X,node1Y,node2X,node2Y)/4;
			}
		}
	}
    return  Be;


}








/**
*Fonction qui assemble la matrice des termes surfaciques de phi*phi de la forme bilineaire
*dans le cas par exemple d une condition de robin fourier 
*/
DenseMatrix assemblageBil(DenseMatrix const& tabNodes, vector<Vecteur> const& tabBound){
   
    
    DenseMatrix B(tabNodes.getnbR(),tabNodes.getnbR());


    for(unsigned int i = 0; i<tabBound.size(); i++){
		
		
		double node1X = tabNodes(tabBound[i](0),0);
		double node1Y = tabNodes(tabBound[i](0),1);
        double node2X = tabNodes(tabBound[i](1),0);
		double node2Y = tabNodes(tabBound[i](1),1);
		
		
       
        DenseMatrix Be = boundElemBil(node1X,node1Y,node2X,node2Y);

		
		for(int j = 0; j<2; j++){
			for(int k = 0; k<2; k++){
				 B(tabBound[i](j),tabBound[i](k)) += Be(j,k);
            }
		}
	
        
	}
	
    return B;
	
	
}

