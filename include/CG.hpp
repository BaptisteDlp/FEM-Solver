#ifndef CG_HPP
#define CG_HPP


#include"Vecteur.hpp"
#include<cmath>
#include<ctime>
#include"mesh.hpp"


template<typename T> class CG; 

template <typename T>
class CG{
	
public:
  CG(T const& m, Vecteur const &v, Vecteur const& g, double tol , int maxit);
		
  void Solve(Mesh const& Th);
  void SolvePrecond(T const& C, Mesh const& Th);
		
  void printDonnees() const;
  void printIter() const;
		
  void NormalisationResidu();
		
  Vecteur getX() const;
  int getnIter() const;
  float getTemps() const;
  double getResVec1(int i) const;
  double getResVec(int i) const;
	
	
private:
  T A;
  Vecteur b;
  Vecteur g0;
  double epsilon;
  int nMax;
  bool resolution;
		
  Vecteur X;
  int nIter;
  std::vector<double> pVec;
  std::vector<double> rVec;
  float tempsExe;
		
	
	
	
	
	
};



template <typename T>
CG<T>::CG(T const& m, Vecteur const& v, Vecteur const& g, double tol , int maxit): A(m), b(v), g0(g), epsilon(tol), nMax(maxit),resolution(false),
											     X(b.getTaille()), nIter(0), tempsExe(0){}
	



	
	





	
template <typename T>
void CG<T>::Solve(Mesh const& Th){
	
  if(!(resolution)){
		
    clock_t t = clock();
		
    Vecteur p(g0);
    Vecteur r(g0);
		
    pVec.push_back(p.norm2());
    rVec.push_back(r.norm2());

		

    double alpha(0);
    double beta(0);
	
		
	       	       	
    while(rVec[nIter] > epsilon*b.norm2() && nIter < nMax){
			
			
      alpha = rVec[nIter]*rVec[nIter]/((A*p)*p);
      X = X + alpha*p;
			
      r = r - alpha*(A*p);
      rVec.push_back(r.norm2());

      for(unsigned int i = 0; i<Th.get_nbvBound(); i++){
	r(Th.get_vertex_bound(i)) = 0;
      }
			
		
				
      beta = rVec[nIter+1]*rVec[nIter+1]/(rVec[nIter]*rVec[nIter]);
      p = r + beta*p;
      pVec.push_back(p.norm2());

			
      nIter++;
			
			
    }
		
	
	
    t = clock() -t;
		
    tempsExe = ((float)t)/CLOCKS_PER_SEC;
		
    resolution = true;
		
  }else{
    std::cout<<"Probleme deja resolu"<<std::endl;
  }
	
	

	
}



template <typename T>
void CG<T>::SolvePrecond(T const& C, Mesh const& Th){
	
  if(!(resolution)){
		
    clock_t t = clock();
		
   
    Vecteur r(g0);
    Vecteur z(C*r);
    Vecteur p(z);
		
    rVec.push_back(r.norm2());
    pVec.push_back(p.norm2());
		

    double alpha(0);
    double beta(0);
    double prodBk(0);
	
    while(rVec[nIter] > epsilon*b.norm2() && nIter < nMax){			
     
      alpha = r*z/(p*(A*p));
      X = X + alpha*p;
			
      prodBk = r*z;
      r = r - alpha*(A*p);
      rVec.push_back(r.norm2());

      for(unsigned int i = 0; i<Th.get_nbvBound(); i++){
	r(Th.get_vertex_bound(i)) = 0;
      }
			
      z = C*r;				
     
      beta = (r*z)/prodBk;
     
      p = z + beta*p;
      pVec.push_back(p.norm2());
			
      nIter++;			
    }
			
    t = clock() -t;
		
    tempsExe = ((float)t)/CLOCKS_PER_SEC;
		
    resolution = true;
		
  }else{
    std::cout<<"Problem already solved"<<std::endl;
  }
	
	

	
}












template <typename T>
void CG<T>::printDonnees() const{
	
  std::cout<<"Resolution of Ax=b with Conjuguate Gradient"<<std::endl;
  std::cout<<"-----Data-----"<<std::endl;
  std::cout<<"A="<<A<<std::endl;
  std::cout<<"b="<<b<<std::endl;
  std::cout<<"tolerance = "<<epsilon<<std::endl;
  std::cout<<"maximum number of iteration = "<<nMax<<"\n"<<std::endl;
	
}






template <typename T>
void CG<T>::printIter() const{
	
  if(resolution){
    std::cout<<"Problem solved with Conjuguate Gradient"<<std::endl;
    if(nIter == nMax){
      std::cout<<"Method did not converge"<<std::endl;
    }else{
      std::cout<<"Method converged with "<<nIter<<" iterations in "<<tempsExe<<"sec"<<std::endl;
    }
  }else{
    std::cout<<"Problem not solved yet"<<std::endl;
  }
	
}







template <typename T>
Vecteur CG<T>::getX() const{ return X; }

template <typename T>  
float CG<T>::getTemps() const{ return tempsExe; }

template <typename T>  
int CG<T>::getnIter() const{ return nIter; } 




#endif
