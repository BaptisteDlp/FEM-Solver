#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<cstring>
#include"Mesh.hpp"







using namespace std;





/**
*Fonction loadNodes
*Fonction qui charge le tableau des noeuds
*/
DenseMatrix loadNodes(string fileName, bool offset){
	
	
  ifstream ifs(fileName);
  double x, y;
  int nbv, label;
	
  ifs>>nbv; ifs>>label;
  if(offset){
    ifs >>label;
  }
  DenseMatrix tabv(nbv,3);
  
  for(int i=0; i<nbv; i++){
    ifs>>x; ifs>>y; ifs>>label;		
    tabv(i,0) = x;
    tabv(i,1) = y;
    tabv(i,2) = label;
  }
  ifs.close();
  return tabv;
				

}




/**
*Fonction loadElements
*Fonction qui charge le tableau des elements
*/
DenseMatrix loadElements(string fileName, bool offset){
	
	
  ifstream ifs(fileName);
  string line;
  double v1, v2, v3, label;
  int nbv, nbt, edge;
       
		
  ifs>>nbv; ifs>>nbt;
  if(offset){
    ifs >> edge;
  }
  
  for(int i=0; i<nbv+1; i++){
    getline(ifs,line); 
  }
		
	
  DenseMatrix tabt(nbt,4);
		
  for(int i=0; i<nbt; i++){
    ifs>>v1; ifs>>v2; ifs>>v3; ifs>>label;			
    tabt(i,0) = v1-1;
    tabt(i,1) = v2-1;
    tabt(i,2) = v3-1;
    tabt(i,3) = label;
  }
		
  ifs.close();
  return tabt;		
	
}



/**
 *Fonction qui creer un fichier de maillage rectangulaire uniforme qui prend en parametrele nombre
 *de noeuds horizontal, vertical, et les les points de depart et d arrivee du domaine
 */
void rectMesh(string fileName, int nbNodesH, int nbNodesV, double x0, double xN, double y0, double yN){
	
	
  ofstream os(fileName);
	
  if(os){
    double nbNodes = nbNodesH*nbNodesV;
    os<<nbNodes<<" ";

    int nbElts = (nbNodesH-1)*(nbNodesV-1)*2;
    os<<nbElts<<endl;

    
    double deltaH = abs(xN-x0)/(nbNodesH-1); double deltaV = abs(yN-y0)/(nbNodesV-1);
		
    vector<double> nodesH; vector<double> nodesV;
		
    for(int i = 0; i<nbNodesH; i++){
      nodesH.push_back(x0+i*deltaH);
    }
    for(int i = 0; i<nbNodesV; i++){
      nodesV.push_back(y0+i*deltaV);
    }
		
    int i = 0;
    int j = 0;
    int k = 1;
	
    while(k<nbNodes+1){
      while(j<nbNodesH){
	os<<nodesH[j]<<" "<<nodesV[i]<<" 0"<<endl;
	j += 1;
	k += 1;
      }
      i += 1;
      j = 0;
    }
	
	
       
    
    int node1 = 1;
    int node2 = 2;
    int node3 = nbNodesH + 2;
    k = 1;
    int tri = 0;

    while(k<nbElts+1){
      tri = 0;
      while(tri<2){
	os<<node1<<" "<<node2<<" "<<node3<<" "<<k<<endl;
					
	if(tri == 0){
	  node2 += nbNodesH;
	  node3 -= 1;
	  k += 1;
	}	
	tri += 1;
      }  
      node1 += 1;
      node2 = node1 + 1;
      node3 += 2;

      if(node1 % nbNodesH == 0){
	node1 += 1;
	node2 += 1;
	node3 += 1;
      }
      k += 1;
                
    }
   			
		
    os.close();
  }else{
    cout<<"ERROR OPENING MESHFILE"<<endl;
  }
	
	
	
}












vector<pair<int,pair<double,double>>>  meshExtraction(string const&  meshNameFile ,string const& regionNameFile, string const& meshLocNameFile, int regionLabelIn, int regionLabelOut, int labelIn, int labelOut, int labelCommon){


  //Opening of the  meshFile
  ifstream os(meshNameFile);
  if(!os){
    cout<<meshNameFile<<" NOT FOUND"<<endl; exit(1);
  }

   
  int nbvG, nbtG;
  int nbvLoc(0), label(0);

  os>>nbvG;
  os>>nbtG;
  os>>label;


 

  //Construction of the array which contains the region label of each vertex
  int tabRegion[nbvG];
  {
    ifstream ifs(regionNameFile);
    if(!ifs){
      cout<<regionNameFile<<" NOT FOUND"<<endl; exit(1);
    }

    for(int i = 0; i<nbvG; i++){
      int region; ifs>>region;
      tabRegion[i] = region;
    }
    ifs.close();
   
  }



  //Opening an auxiliary local meshfile
  string mesh("mesh.msh");
  ofstream ofs(mesh);
  if(!ofs){
    cout<<mesh<<" NOT FOUND"<<endl; exit(1);
  }

  

    
  //Initialization of the array which states whether a vertex is part of the new numerotation  or not
  int  tabInd[nbvG];
  for(int i = 0; i < nbvG; i++){
    tabInd[i] = -1;
  }

  vector<pair<int,pair<double,double>>> gamma;
 

  //Construction of the new vertices numerotation
  for(int i = 0; i<nbvG; i++){
    double x,y;
    os>>x; os>>y; os>>label;
    if(tabRegion[i] == regionLabelIn || tabRegion[i] == regionLabelOut || label == labelIn || label == labelOut || label == labelCommon){
      tabInd[i] = ++nbvLoc;
      ofs<<x<<" "<<y<<" "<<" "<<label<<endl;
      if(label == labelIn || label == labelCommon){
	//add this vertex to gamma
	pair<double,double> coord(x,y);
	pair<int,pair<double,double>> p(nbvLoc-1,coord);
	gamma.push_back(p);
      }
    }
  }

  

  
  int v1, v2, v3;
  int nbtLoc(0);


  //Construction of the new elements numerotation
  for(int i = 0; i<nbtG; i++){
    os>>v1; os>>v2; os>>v3; os>> label;
    if(tabInd[v1-1] != -1 && tabInd[v2-1] != -1 && tabInd[v3-1] != -1){
      nbtLoc++;
      ofs<<tabInd[v1-1]<<" "<< tabInd[v2-1]<<" "<<tabInd[v3-1]<<" "<< label<<endl;
    }
  }


  
 
  ofs.close();
  os.close();
  cout<<nbvLoc<<" "<<nbtLoc<<endl;
  

  //Creation of the final local mesh file
  {
    ofstream file("file.txt");
    file<<nbvLoc<<" "<<nbtLoc<<endl;
    string s;
    s = "cat file.txt " + mesh + " > " + meshLocNameFile;
    char * cstr = new char [s.length()+1];
    strcpy (cstr, s.c_str());
    system(cstr);
    delete [] cstr;
    file.close();
    system("rm file.txt");
    system("rm mesh.msh");
  }



  return gamma;


}






vector<pair<int,int>> buildGammaij(DenseMatrix const& tabv1, DenseMatrix const& tabv2,int labelIn){

  vector<pair<int,int>> gamma;
  
  for(int i = 0; i<tabv1.getnbR(); i++){
    if(tabv1(i,2) == labelIn){
      for(int j = 0; j<tabv2.getnbR(); j++){
	if(tabv2(j,2) == labelIn){
	  if(tabv2(j,0) == tabv1(i,0) && tabv2(j,1) == tabv1(i,1)){
	    pair<int,int> p(i,j);
	    gamma.push_back(p);  break;
	  }
	}

      }
    }
  }

  return gamma;
}








vector<pair<int,int>> buildGamma(vector<pair<int,pair<double,double>>> const& g, DenseMatrix const& tabv1,  DenseMatrix const& tabv2,int labelIn){

  vector<pair<int,int>> gamma;

  for(unsigned int i = 0; i<g.size(); i++){
    int ind = g[i].first;
    double x = g[i].second.first;
    double y = g[i].second.second;
    for(int j = 0; j<tabv2.getnbR(); j++){
      if(tabv2(j,2) == labelIn){
	if(tabv2(j,0) == tabv1(ind,0) && tabv2(j,1) == tabv1(ind,1)){
	  // if(x == tabv2(j,0)  && y == tabv2(j,1)){ //this test does not work, I do not know why ????
	  pair<int,int> p(ind,j);
	  //cout<<tabv2(j,0)<<" "<<"=="<< g[i].second.first<<" && "<<tabv2(j,1)<<"=="<< g[i].second.second<<endl;
	  gamma.push_back(p);  break;
	}
      }

    }
    
  }
 
  return gamma;

}









void exportMeshData(string meshName, DenseMatrix const& tabNodes,DenseMatrix const& tabElts){
  
  ofstream os(meshName);

  for(int i = 0; i<tabElts.getnbR(); i++){
    for(int j = 0; j<3; j++){
      os<<tabNodes(tabElts(i,j),0)<<" "<<tabNodes(tabElts(i,j),1)<<" "<<tabNodes(tabElts(i,j),2)<<endl;
    }
    os<<tabNodes(tabElts(i,0),0)<<" "<<tabNodes(tabElts(i,0),1)<<" "<<tabNodes(tabElts(i,0),2)<<"\n"<<endl;
  }
  os.close();

}








	
void exportMeshGnuplot(string scriptName,string meshName){
  
  ofstream os(scriptName);
	
  os<<"set title \"Mesh\""<<endl;
  os<<"set key off"<<endl;
  os<<"plot '"<<meshName<<"' with lines"<<endl;
  os.close();

       	
}



	
	
void exportSolData(string SolName, DenseMatrix const& tabNodes, DenseMatrix const& tabElts, Vecteur const& solution){
	
  ofstream os(SolName);
	
  if(os){
    for(int i = 0; i<tabElts.getnbR(); i++){
      for(int j = 0; j<3; j++){
	os<<tabNodes(tabElts(i,j),0)<<" "<<tabNodes(tabElts(i,j),1)<<" "<<solution(tabElts(i,j))<<endl;
	//os<<"\n";
      }
      os<<tabNodes(tabElts(i,0),0)<<" "<<tabNodes(tabElts(i,0),1)<<" "<<solution(tabElts(i,0))<<"\n\n"<<endl;
    }
    os.close();
  }else{
    cout<<"ERROR OPENING DATASOLUTION"<<endl;
  }
	
}





void exportSolGnuplot(string scriptName,string SolName, string title){
	
  ofstream os(scriptName);
	
  if(os){
    os<<"set title \""<<title<<"\""<<endl;
    os<<"set key off"<<endl;
    os<<"splot'"<<SolName<<"' using 1:2:3 with lines palette"<<endl;
    os.close();
  }else{
    cout<<"ERROR OPENING SOLUTIONSCRIPT"<<endl;
  }
	
}






  












vector<vector<double>> sol(DenseMatrix const& tabNodes, Vecteur const& solution){
	
  vector<vector<double>> tab;
	
  for(int i = 0; i<tabNodes.getnbR(); i++){
    vector<double> v;
    v.push_back(tabNodes(i,0)); v.push_back(tabNodes(i,1)); v.push_back(solution(i));
    tab.push_back(v);
  }
	
  return tab;
	
}










void exportSolGnuplotPm3d(std::string scriptName, std::string SolName, std::string title){
	
	
  ofstream os(scriptName);
	
  if(os){
    os<<"set title \""<<title<<"\""<<endl;
    os<<"set pm3d at ss"<<endl;
    os<<"set hidden3d"<<endl;
    os<<"set surface"<<endl;
    os<<"set contour base"<<endl;
    os<<"set view 0,0"<<endl;
    os<<"set key off"<<endl;
    os<<"splot 'C:/Users/Delaporte/Desktop/FEMSolver/"<<SolName<<"' using 1:2:3 with lines"<<endl;
    os.close();
  }else{
    cout<<"ERROR OPENING SOLUTIONSCRIPT PM3D"<<endl;
  }
	
	
	
	
	
}

