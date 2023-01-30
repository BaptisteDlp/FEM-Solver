#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<cstring>
#include<cmath>
#include"mesh.hpp"


using namespace std;


void  build_edges(double **tabv, int nbv, int nbt, int *nbe, int *nbeb, int *nbvBound, int **tabt, vector<std::pair<int,int>> &tabe, vector<std::pair<int,int>> &tabeb, vector<int> &tabvBound){

  int *First = new int[nbv];
  for(int i = 0; i < nbv; i++){
    First[i] = -1;
  }

  int *element = new int[nbt*3]; //element[k+j] == edge j of element k
  int *element_adj = new int [nbt*3]; //element_adj[k+j] == edge j of adjacent element of k
  
  for(int i = 0; i < nbt*3; i++){
    element_adj[i] = -1;
  }

  int *num_P2 = new int [nbt*6];
  
  for(int i = 0; i < nbt*6; i++){
    num_P2[i] = -3;
  }

  int ndof_P2 = 0;

  vector<int> Next;
  vector<int> occ;

  
  //loop over every elements
  for(int t = 0; t < nbt; ++t){
    //loop over eveery edges of element t
    cout<<"triangle "<< t<<endl;
    for(int j = 0; j<3; ++j){

      int  k1 = (j+1) % 3;
      int  k2 = (j+2) % 3;
      int il = min(tabt[t][k1],tabt[t][k2]);
      int jl = max(tabt[t][k1],tabt[t][k2]);

      //======DEBUG======
      cout<<"arete "<<j<<endl;
      cout<<"k1="<<k1<<" k2="<<k2<<endl;
      cout<<"il="<<il<<" jl="<<jl<<endl;
      cout<<"("<<tabv[il][0]<<","<<tabv[il][1]<<")-----("<<tabv[jl][0]<<","<<tabv[jl][1]<<")"<<endl;

  
      bool exist(false);

      int k = First[il];
    
      while(k != -1 && !(exist)){

	if(tabe[k].first == il && tabe[k].second == jl){
	  exist = true; 
	}

	if(!(exist)){
	  k = Next[k];
	}

      }
     
      if(!(exist)){
	//edge  j does not exist, we have to add it
	pair<int,int> edge;
	edge.first = il;
	edge.second = jl;
	tabe.push_back(edge);

	element[*nbe] = 3*t + j;
	
	
	Next.push_back(First[il]);
	First[il] = *nbe;
	occ.push_back(1);
	(*nbe)++;
      }else{
	//edge j already exists
	occ[k] ++;

	int kk = element[k]/3, aa = element[k]%3;
	//========DEBUG========
	cout<<"t="<<t<<"|j="<<j<<endl;
	cout<<"kk="<<kk<<"|aa="<<aa<<endl;
	element_adj[kk*3+aa] = t*3+j; //update adjacent element
	element_adj[t*3+j] = kk*3+aa;
	element[t*3+j] = kk*3+aa;
	

	
      }   
    }
    cout<<endl;//DEBUG===================
  }




  
  //======================DEBUG=====================//
  cout<<"element"<<endl;
  for(int i = 0; i<nbt; i++){
    for(int j = 0; j < 3; j++){
      cout<<element[3*i+j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;

  cout<<"element_adj"<<endl;
  for(int i = 0; i<nbt; i++){
    for(int j = 0; j < 3; j++){
      cout<<element_adj[3*i+j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
  //==========================END DEBUG=====================//





  
  int *dof_vertex = new int[nbv]; // dof_vertex[s] => number of vertex s degree of freedom
  for(int i = 0; i < nbv; i++){
    dof_vertex[i] = -1;
  }
  int knum =0;

  //numerotation P2
    
  //loop over every elements
  for(int t = 0; t < nbt; t++){
    //loop over every vertices
    for(int i=0; i<3; ++i){
      int s = tabt[t][i]; 
      if(dof_vertex[s] < 0){ // dof of vertex s does not exist, then  dof_vertex[s]<-ndof;
	dof_vertex[s] = ndof_P2++;
      }
					
      num_P2[knum++] = dof_vertex[s]; //update numerotation
    }
    knum = knum + 3;
  }
    
  /*A ce stade, on a construit les ddl de chaque sommet de chaque triangle t*/
  /*Reste a construire ceux associes aux aretes*/
			
  //loop over every edges of element t
  knum = 3;
  for(int t = 0; t < nbt; t++){
    for(int  j = 0; j < 3; ++j){
      int kk = element_adj[3*t+j] / 3; //On recupere le numero du triangle adjacent a k
      int aa = element_adj[3*t+j] % 3; //On recupere le numero de l arete commune 
				
      if(element_adj[3*t+j] >=0  && kk < t){ // on a deja construit le numero
	num_P2[knum++] = num_P2[6*kk+3+aa]; //offset =  6*nb element + 3*vertex + aa
      }else{
	num_P2[knum++] = ndof_P2++;
      }
    }
    knum = knum + 3;
  }
  
  assert(ndof_P2 == *nbe + nbv);

  for(int k=0; k<nbt; ++k){
    cout << k <<" : ";
    for(int i=0; i<6; ++i)
      cout << num_P2[6*k+i]<< ' ';
    cout << endl;
				 
  }
    
  
 for(int q = 0; q<*nbe; q++){
    if(occ[q] == 1){
      pair<int,int> edge;
      edge.first = tabe[q].first;
      edge.second = tabe[q].second;
      tabeb.push_back(edge); (*nbeb) ++;	
    }
  }
  cout<<"\t-->edges loaded"<<endl;


  for(int i = 0; i < nbv; i++){
    First[i] = -1;
  }
   

  //loop over every edges on the boundary
  for(int i = 0; i < tabeb.size(); i++){
    
    int v1 = tabeb[i].first;
    int v2 = tabeb[i].second;

    if(First[v1] == -1){
      //vertex v1 does not exist, we have to add it
      tabvBound.push_back(v1);
      First[v1] =  1; (*nbvBound)++;
    }
		
    if(First[v2] == -1){
      //vertex v2 does not exist, we have to add it
      tabvBound.push_back(v2);
      First[v2] =  1; (*nbvBound)++;
    }
  }
  cout<<"\t-->vertices on boundary loaded"<<endl;

  delete [] First;
  delete [] element;
  delete [] element_adj;
  delete [] dof_vertex;
  delete [] num_P2;
}




Mesh::Mesh(std::string file_name, bool offset){


  cout<<"==Mesh loading=="<<endl;
  mesh_name = file_name;
  ifstream file(file_name);
  if(!file){
    cout<<"ERROR OPENNING "<<file_name<<endl; exit(1);
  }

  int var1(0);
  double var2(0);
  file >> nbv; file >> nbt;
  if(offset){
    file >> nbe;
  }

  tabv = new double*[nbv];
 
  //loading vertices
  for(int i = 0; i < nbv; i++){
    tabv[i] = new double[2];
    file >> tabv[i][0];
    file >> tabv[i][1];
    file >> var1;
  }
  cout<<"\t-->vertices loaded"<<endl;
  

  tabt = new int*[nbt];

  //loading elements
  for(int i = 0; i < nbt; i++){
    tabt[i] = new int[4];
    file >> tabt[i][0]; tabt[i][0]--; //vertex 1
    file >> tabt[i][1]; tabt[i][1]--; //vertex 2
    file >> tabt[i][2]; tabt[i][2]--; //vertex 3
    file >> tabt[i][3]; //label
  }
  cout<<"\t-->elements loaded"<<endl;
 
  //loading edges 
  
  nbe = 0; nbeb = 0; nbvBound = 0;
  build_edges(tabv,nbv,nbt,&nbe,&nbeb,&nbvBound,tabt,tabe,tabeb,tabvBound);

 

  

  mesh_data = "data/mesh.data";
  mesh_gnu = "data/mesh.gnu";
  solution_data = "data/solution.data";
  solution_gnu = "data/solution.gnu";
  resolved = false;
  
  file.close();

  cout<<"\t"<<nbv<<" vertices"<<endl;
  cout<<"\t"<<nbt<<" elements"<<endl;
  cout<<"\t"<<nbe<<" edges"<<endl;
  cout<<"\t"<<nbvBound<<" vertices on boundary"<<endl;
  cout<<"==Mesh loaded=="<<endl;
  


}





Mesh::~Mesh(){

  for(int i = 0; i < nbv; i++){
    delete [] tabv[i];
  }

  delete [] tabv;

  for(int i = 0; i < nbt; i++){
    delete [] tabt[i];
  }

  delete [] tabt;

}





void Mesh::uniform_mesh(int nbv_h, int nbv_v, double x0, double xN, double y0, double yN, std::string fileName){
	
	
  ofstream os(fileName);
  if(!os){
    cout<<"ERROR OPENNING "<<fileName<<endl; exit(1);
  }
	

  double nbv = nbv_h*nbv_v;
  int nbt = (nbv_h-1)*(nbv_v-1)*2;
  os<<nbv<<" "<<nbt<<endl;

    
  double deltaH = abs(xN-x0)/(nbv_h-1);
  double deltaV = abs(yN-y0)/(nbv_v-1);
		
  vector<double> nodesH; vector<double> nodesV;
		
  for(int i = 0; i<nbv_h; i++){
    nodesH.push_back(x0+i*deltaH);
  }
  for(int i = 0; i<nbv_v; i++){
    nodesV.push_back(y0+i*deltaV);
  }
		
  int i(0), j(0), k(1);
	
  while(k<nbv+1){
    while(j<nbv_h){
      os<<nodesH[j]<<" "<<nodesV[i]<<" 0"<<endl;
      j ++;
      k ++;
    }
    i ++;
    j = 0;
  }
	
  int v1(1), v2(2), v3(nbv_h+2), tri(0);
  k = 1;

  while(k<nbt+1){
    tri = 0;
    while(tri<2){
      os<<v1<<" "<<v2<<" "<<v3<<" "<<k<<endl;
					
      if(tri == 0){
	v2 += nbv_h;
	v3 --;
	k ++;
      }	
      tri ++;
    }  
    v1 ++;
    v2 = v1 + 1;
    v3 += 2;

    if(v1 % nbv_h == 0){
      v1 ++;
      v2 ++;
      v3 ++;
    }
    k ++;
                
  }
   			
		
  os.close();
 
	
	
	
}















void Mesh::print_vertices() const {

  cout<<"\ntabv"<<endl;
  for(int i = 0; i < nbv; i++){
    cout<<"("<<tabv[i][0]<<","<<tabv[i][1]<<") "<<endl;
  }

}

void Mesh::print_elements() const {

  cout<<"\ntabt"<<endl;
  for(int i = 0; i < nbt; i++){
    cout<<"|"<<tabt[i][0]<<","<<tabt[i][1]<<","<<tabt[i][2]<<"| "<<tabt[i][3]<<endl;
  }
}


void Mesh::print_edges() const {

  cout<<"\ntabe"<<endl;
  for(int i = 0; i < nbe; i++){
    cout<<"|"<<tabe[i].first<<","<<tabe[i].second<<"|"<<endl;
  }

}

void Mesh::print_vertices_boundary() const {

  cout<<"\ntabvBound"<<endl;
  for(int i = 0; i < nbvBound; i++){
    cout<<"|"<<tabvBound[i]<<"|"<<endl;
  }

}



void Mesh::export_mesh_data() const{

  {
    ofstream file(mesh_data);

    for(int i = 0; i<nbt; i++){
      for(int j = 0; j<3; j++){
	file << tabv[tabt[i][j]][0] << " " << tabv[tabt[i][j]][1] << " 0"<<endl;
      }
      file << tabv[tabt[i][0]][0] << " " << tabv[tabt[i][0]][1] << " 0\n" <<endl;
    }
  
    file.close();
  }

  {
    ofstream file(mesh_gnu);

    file<<"set title \"Mesh\""<<endl;
    file<<"set key off"<<endl;
    file<<"plot '"<<mesh_data<<"' with lines"<<endl;
    file.close();
  
    file.close();
  }

}


void Mesh::mesh_visu(bool mesh, bool solution) const{
  string name;
  if(mesh){
    name = mesh_gnu;
  }

  if(solution){
    name = solution_gnu;
  }
  
  string s = "gnuplot  " + name + " -persist";
  char * cstr = new char [s.length()+1];
  strcpy (cstr, s.c_str());
  system(cstr);
 
}



void Mesh::export_solution_data(Vecteur const& solution) const{

  if(!resolved){
    cout<<"PROBLEM NOT SOLVED YET"<<endl; exit(1);
  }

  {
    ofstream file(solution_data);

    for(int i = 0; i<nbt; i++){
      for(int j = 0; j<3; j++){
	file << tabv[tabt[i][j]][0] << " " << tabv[tabt[i][j]][1] << " "<<solution(tabt[i][j])<<endl;
      }
      file << tabv[tabt[i][0]][0] << " " << tabv[tabt[i][0]][1] << " "<<solution(tabt[i][0])<<endl;
    }
  
    file.close();
  }

  {
    ofstream file(solution_gnu);

    file<<"set title \"Approximate solution\""<<endl;
    file<<"set key off"<<endl;
    file<<"splot '"<<solution_data<<"'using 1:2:3 with lines palette lw 2"<<endl;
    file.close();
  
    file.close();
  }

}





int Mesh::get_nbv() const { return nbv;} 
int Mesh::get_nbt() const { return nbt;} 
int Mesh:: get_nbe() const {return nbe;} 
int Mesh::get_nbeb() const {return nbeb;} 
int Mesh::get_nbvBound() const {return nbvBound;} 
double Mesh::get_x(int index) const {assert(index >=0 && index < nbv); return tabv[index][0];} 
double Mesh::get_y(int index) const {assert(index >=0 && index < nbv); return tabv[index][1];}
int Mesh::get_vertex(int index, int v) const {assert(index >=0 && index < nbt && v >= 0 && v < 3); return tabt[index][v];} 
int Mesh::get_vertex_bound(int index) const {assert(index >=0 && index < nbv); return tabvBound[index];}

/*setters*/
void Mesh::set_resolved(){ resolved = true; }
void Mesh::unset_resolved(){ resolved = false; }

int Mesh::get_edge(int index, int v) const{
  assert(index >=0 && index < nbe && v >= 0 && v < 2);
  double res;
  if(v == 0){
    res = tabe[index].first;
  }else{
    res = tabe[index].second;
  }
  return res;
}



