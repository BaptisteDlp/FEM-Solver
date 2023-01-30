#include<iostream>
#include<cstdlib>
#include<cassert>
#include<cmath>
#include<ctime>
#include<fstream>
#include"DenseMatrix.hpp"
#include"SparseMatrixCoo.hpp"
#include"Assembly.hpp"
#include"Penalization.hpp"
#include"Vecteur.hpp"
#include"Mesh.hpp"
#include"mesh.hpp"
#include"Edge.hpp"
#include"CG.hpp"








using namespace std;




Vecteur operator*(MatCSR const& M, Vecteur const& v){
	
  assert(M.m == v.getTaille());
	
  Vecteur Mv(v.getTaille());
	
  return M.addmatmul(v,Mv);
	
}



void printMatCSR1(MatCSR const& A){
  for(int i=0; i<A.n; ++i)
    {
      std::cout<< i << " : " ;
        for(int k=A.ip[i]; k <A.ip[i+1]; ++k)
	  std::cout<< " " << A.j[k] << " " << A.a[k] << ";   ";
	std::cout << std::endl;
    }

}







void pseudoElimination(MatCSR & A, Mesh const& Th){


  for(unsigned int i = 0; i<Th.get_nbvBound(); i++){
     

    int il = Th.get_vertex_bound(i);

    for(int k = A.ip[il]; k<A.ip[il+1]; ++k){
      if(A.j[k] == il){
	 A(il,A.j[k]) = 1;
      }else{
	A(il,A.j[k]) = 0;
       	A(A.j[k],il) = 0;
      }
    }
             
    
  }


  
}



Vecteur setFunction(int ndof, Mesh const& Th){
  Vecteur f(ndof);
  for(int i = 0; i<ndof; i++){
    double x = Th.get_x(i);
    double y = Th.get_y(i);
    f(i) = 1;
  }
  return f;
}




Vecteur setbi(MatCSR const& M, Vecteur const& f, Mesh const& Th){
  Vecteur bi = M*f;
  for(unsigned int i = 0; i<Th.get_nbvBound(); i++){
    bi(Th.get_vertex_bound(i)) = 0;
  }
  
  return bi;
}










int main(int argc, char **argv){
	
  if(argc<=1){
    cout<<"Please precise mesh file"<<endl;
    cout<<"./bin/poisson_solver data/Th.msh"<<endl;exit(1);
  }
  

 
  cout<<"################"<<endl;
  cout<<"#Poisson solver#"<<endl;
  cout<<"################\n"<<endl;
	
	
  string Omega(argv[1]);
 
    
  Mesh Th(Omega);

  Th.uniform_mesh(2,2,0,0.5,0,0.5);

  Th.export_mesh_data();
  Th.mesh_visu();

  Th.print_vertices();
  Th.print_elements();
  Th.print_edges();

  
  int ndof = Th.get_nbv();
 

  /*Pre-processing assembly*/

  clock_t t = clock();
  MatCSR M = Stiffness_Mass_matrix(Th,false); //Mass Matrix
  MatCSR K = Stiffness_Mass_matrix(Th); //Stiffness matrix
  pseudoElimination(K,Th);
  Vecteur f = setFunction(ndof,Th);//set 1 to every entries of f
  Vecteur b = setbi(M,f,Th);//set the initial second member of the linear system: M*f and 0 for the entries which correspond to vertices on boundary
  /*build Jacobi preconditionner for CG*/
  MatrixMap C_map;

  for(int i = 0; i<ndof;i++){
    pair<int,int> pos(i,i);
    C_map[pos] = 1/K(i,i);
  }
  MatCSR C(ndof,ndof,C_map);
  Vecteur g0 = b; //Initiazation of first residual part for CG
  t = clock() -t;
  
  float execTimeA = ((float)t)/CLOCKS_PER_SEC;

   	
  t = clock();
  cout<<"\n=====LINEAR SYSTEM SOLVING====="<<endl;
  double tol(0.001);
  int maxIt(20000);
  CG<MatCSR> System(K,b,g0,tol,maxIt);
  System.SolvePrecond(C,Th);//Solve with preconditionner C
  System.printIter();
  Vecteur solution = System.getX();
  t = clock() -t;


  Th.set_resolved();
  Th.export_solution_data(solution);
  Th.mesh_visu(false,true);

  
  float execTimeR = ((float)t)/CLOCKS_PER_SEC;
  cout<<"\nCPU TIME for Assembly:"<<execTimeA<<"sec"<<endl;
  cout<<"\nCPU TIME for resolution:"<<execTimeR<<"sec"<<endl;

   
  double maxV(0);
 
  for(int i = 0; i<solution.getTaille(); i++){
    maxV = max(maxV,solution(i));
  }
 
  cout<<"max=" <<maxV<<endl;

   
  return 0;
}
