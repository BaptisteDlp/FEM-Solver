#include<iostream>
#include"Edge.hpp"
#include<vector>

using namespace std;



/**
 *Fonction qui renvoie le tableau des aretes
 */
vector<Vecteur> findEdge(DenseMatrix const& tabNodes, DenseMatrix const& tabElts){
    
  vector<Vecteur> Edge;
  Vecteur First(tabNodes.getnbR());
  First.ones();
  First = (-1)*First;
	
	
  vector<int> Next;
  int edge1, edge2, aux;

  int nbEdges = 0;
  int i = 0;
                         
  while(i<tabElts.getnbR()){
    //on parcourt tout les triangles
        
    for(int j = 0; j<3; j++){

      //on parcourt toutes les aretes du triangle
      if(j == 2){
	edge1 = tabElts(i,0);
	edge2 = tabElts(i,2);
      }else{
	edge1 = tabElts(i,j);
	edge2 = tabElts(i,j+1);
      }


      if(edge2 < edge1){
	//on rearrange l ordre des noeuds (croissant)
	aux = edge2;
	edge2 = edge1;
	edge1 = aux;
      }
                
           
            
      if(First(edge1) == -1){
	//l arete n existe pas, on la rajoute
                
	Vecteur edge(2);
	edge(0) = edge1; edge(1) = edge2;
	Edge.push_back(edge);
                
				
	//mise a jour des autres tableaux
	Next.push_back(First(edge1));
	First(edge1) = nbEdges;
	nbEdges += 1;
      }else{

	//il faut verifier si l arete existe deja ou pas
	int k = First(edge1);
                
	int exist = 0;
               
	while(k != -1 && exist ==0){
	  //on parcourt la classe d equivalence
	  if(Edge[k](0) == edge1 and Edge[k](1) == edge2){
	    exist = 1;
	  }
                          
	  if(exist == 0){
	    k = Next[k];
	  }
	}
                         
	if(exist == 0){
	  //l arrete n existe pas, on la rajoute
	  Vecteur edge(2);
	  edge(0) = edge1; edge(1) = edge2;
	  Edge.push_back(edge);
                
				
	  //#mise a jour des autres tableaux
	  Next.push_back(First(edge1));
	  First(edge1) = nbEdges;
	  nbEdges += 1;
	}
      }
              
    }			  
    i += 1;
  }
	
  return Edge;
	
}







/**
 *Fonction qui renvoie le tableau des aretes aux bords
 */
vector<Vecteur> findBoundary(DenseMatrix const& tabNodes, DenseMatrix const& tabElts){
    
  vector<Vecteur> Edge;
  vector<Vecteur> Boundary;
  Vecteur First(tabNodes.getnbR());
  First.ones();
  First = (-1)*First;
	
  vector<int> occ;
  vector<int> Next;
  int edge1, edge2, aux;

  int nbEdges = 0;
  int i = 0;
                         
  while(i<tabElts.getnbR()){
    //on parcourt tout les triangles
        
    for(int j = 0; j<3; j++){

      //on parcourt toutes les aretes du triangle
      if(j == 2){
	edge1 = tabElts(i,0);
	edge2 = tabElts(i,2);
      }else{
	edge1 = tabElts(i,j);
	edge2 = tabElts(i,j+1);
      }


      if(edge2 < edge1){
	//on rearrange l ordre des noeuds (croissant)
	aux = edge2;
	edge2 = edge1;
	edge1 = aux;
      }
                
           
            
      if(First(edge1) == -1){
	//l arete n existe pas, on la rajoute
                
	Vecteur edge(2);
	edge(0) = edge1; edge(1) = edge2;
	Edge.push_back(edge);
                
				
	//mise a jour des autres tableaux
	Next.push_back(First(edge1));
	First(edge1) = nbEdges;
	occ.push_back(1);
	nbEdges += 1;
      }else{

	//il faut verifier si l arete existe deja ou pas
	int k = First(edge1);
                
	int exist = 0;
               
	while(k != -1 && exist ==0){
	  //on parcourt la classe d equivalence
	  if(Edge[k](0) == edge1 and Edge[k](1) == edge2){
	    exist = 1;
	  }
                          
	  if(exist == 0){
	    k = Next[k];
	  }
	}
                         
	if(exist == 0){
	  //l arrete n existe pas, on la rajoute
	  Vecteur edge(2);
	  edge(0) = edge1; edge(1) = edge2;
	  Edge.push_back(edge);
                
				
	  //#mise a jour des autres tableaux
	  Next.push_back(First(edge1));
	  First(edge1) = nbEdges;
	  occ.push_back(1);
	  nbEdges += 1;
	}else{
	  occ[k] += 1;
	}
				
      }
              
    }			  
    i += 1;
  }
	
	
	
  for(unsigned int q = 0; q<Edge.size(); q++){
    if(occ[q] == 1){
      Vecteur bound(2);
      bound(0) = Edge[q](0); bound(1) = Edge[q](1);
      Boundary.push_back(bound);
    }
  }
	
	
  return Boundary;
	
   
	
}



/*
 *Fonction qui renvoie le tableau des noeuds aux bords
 */
vector<int> findNodesBoundary(DenseMatrix const& tabNodes, vector<Vecteur> const& tabBoundary){

  vector<int> NodesBoundary;
  Vecteur First(tabNodes.getnbR());
  First.ones();
  First = (-1)*First;
   
 
  unsigned int i = 0;
                         
  while(i<tabBoundary.size()){
    //on parcourt toutes les aretes
    int node1 = tabBoundary[i](0);
    int node2 = tabBoundary[i](1);

    if(First(node1) == -1){
      //le noeud n existe pas, on le rajoute
      NodesBoundary.push_back(node1);
      First(node1) =  1;
    }
		
    if(First(node2) == -1){
      //le noeud n existe pas, on le rajoute
      NodesBoundary.push_back(node2);
      First(node2) =  1;
    }
		
    i += 1;
  }
	
  return NodesBoundary;

}







vector<int> idBord(DenseMatrix const& tabNodes, vector<int> const& tabNodesBound){

  //On initialise le tableau a -1
  vector<int> tab(tabNodes.getnbR(),-1);
	
	

  //On parcourt le tableau des noeuds du bords
  for(unsigned int i = 0; i<tabNodesBound.size(); i++){
         
    int k = tabNodesBound[i];
	
       
    if(tabNodes(k,0) == 1){
      //Dans ce cas le noeud k appartient a Gamma1
      int trouve = 0;
      int j = 0;
			
      while(trouve == 0){
	//On cherche le noeud complementaire dans Gamma0
	//parmi les noeuds du bords
	int l = tabNodesBound[j];		 
	if(tabNodes(l,0) == 0 && tabNodes(k,1) == tabNodes(l,1)){
	  //l est le noeud recherche
	  trouve = 1;
	  tab[k] = l;
	}
	j += 1;
      }
			
    }
				 
  }
  
  return tab;

}









 
// Fusion des listes t(de1..vers1) et t(de2..vers2) dans tmp(posInTmp..posInTmp+count-1)












/*
 *Fonction qui affiche le tableau des aretes
 */
void printEdges(vector<Vecteur> const& edge){
	
  for(unsigned int i = 0; i<edge.size(); i++){
    cout<<"| "<<edge[i](0)<<"|"<<edge[i](1)<<" |"<<endl;
  }
	
}

/*
 *Fonction qui affiche le tableau des aretes du bord
 */
void printBound(vector<Vecteur> const& bound){
	
  for(unsigned int i = 0; i<bound.size(); i++){
    cout<<"| "<<bound[i](0)<<"|"<<bound[i](1)<<" |"<<endl;
  }
	
}


/*
 *Fonction qui affiche le tableau des noeuds du bord
 */
void printNodesBound(vector<int> const& nodesBound){
	
  for(unsigned int i = 0; i<nodesBound.size(); i++){
    cout<<"| "<<nodesBound[i]<<" |"<<endl;
  }
	
}
	

	
/*
 *Fonction qui affiche le tableau idBord
 */
void printIdBord(std::vector<int> const& idBord){
	
  for(unsigned int i = 0; i<idBord.size(); i++){
    cout<<"| "<<idBord[i]<<" |"<<endl;
  }
	
}
	
