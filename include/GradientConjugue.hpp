#ifndef GRADIENTCONJUGUE_HPP
#define GRADIENTCONJUGUE_HPP

#include<SparseMatrixCoo.hpp>
#include<Vecteur.hpp>
#include<cmath>


template<typename T> class GradientConjugue; 

template <typename T>
class GradientConjugue{
	
	public:
		GradientConjugue(T const& m, Vecteur const &v, double tol , int maxit);
		
		void Solve();
		void SolvePrecond(T const& P);
		
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
		double epsilon;
		int nMax;
		bool resolution;
		
		Vecteur X;
		int nIter;
		std::vector<double> resVec1;
		std::vector<double> resVec;
		float tempsExe;
		
	
	
	
	
	
};



template <typename T>
GradientConjugue<T>::GradientConjugue(T const& m, Vecteur const& v, double tol , int maxit): A(m), b(v), epsilon(tol), nMax(maxit),resolution(false),
				X(b.getTaille()), nIter(0), tempsExe(0){}
	



	
	





	
template <typename T>
void GradientConjugue<T>::Solve(){
	
	if(!(resolution)){
		
		clock_t t = clock();
		
		Vecteur r1(b);
		Vecteur r(b);
		
		resVec1.push_back(r1.norm2());
		resVec.push_back(r.norm2());

		

		double alpha(0);
		double beta(0);
	
		
		
		
	
		while(resVec[nIter] > epsilon*b.norm2() && nIter < nMax){
			
			//std::cout<<resVec2[nIter]<<std::endl;
			alpha = resVec[nIter]*resVec[nIter]/((A*r1)*r1);
			X = X + alpha*r1;
			
			r = r - alpha*(A*r1);
			resVec.push_back(r.norm2());
			
		
				
			beta = resVec[nIter+1]*resVec[nIter+1]/(resVec[nIter]*resVec[nIter]);
			r1 = r + beta*r1;
			resVec1.push_back(r1.norm2());

			
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
void GradientConjugue<T>::SolvePrecond(T const& P){
	
	if(!(resolution)){
		
		clock_t t = clock();
		
		//Vecteur Pb = P*b;
		//Vecteur r1(Pb);
		//Vecteur r(Pb);
		
		Vecteur r(b);
		Vecteur z(P*r);
		Vecteur r1(z);
		
		resVec1.push_back(r1.norm2());
		resVec.push_back(r.norm2());

		

		double alpha(0);
		double beta(0);
		double prodBk(0);
	
		
		
		
	
		while(resVec[nIter] > epsilon*b.norm2() && nIter < nMax){
			
			//std::cout<<resVec2[nIter]<<std::endl;
			//alpha = resVec[nIter]*resVec[nIter]/(r1*(P*(A*r1)));
			alpha = r*z/(r1*(A*r1));
			X = X + alpha*r1;
			
			//r = r - alpha*(P*(A*r1));
			prodBk = r*z;
			r = r - alpha*(A*r1);
			resVec.push_back(r.norm2());
			
			z = P*r;
			
		
				
			//beta = resVec[nIter+1]*resVec[nIter+1]/(resVec[nIter]*resVec[nIter]);
			beta = (r*z)/prodBk;
			//r1 = r + beta*r1;
			r1 = z + beta*r1;
			resVec1.push_back(r1.norm2());

			
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
void GradientConjugue<T>::printDonnees() const{
	
	std::cout<<"Resolution de Ax=b par GradientConjugue"<<std::endl;
	std::cout<<"-----Donnees du probleme-----"<<std::endl;
	std::cout<<"A="<<A<<std::endl;
	std::cout<<"b="<<b<<std::endl;
	std::cout<<"tolerance = "<<epsilon<<std::endl;
	std::cout<<"nombre iteration max = "<<nMax<<"\n"<<std::endl;
	
}






template <typename T>
void GradientConjugue<T>::printIter() const{
	
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
Vecteur GradientConjugue<T>::getX() const{ return X; }

template <typename T>  
float GradientConjugue<T>::getTemps() const{ return tempsExe; }

template <typename T>  
int GradientConjugue<T>::getnIter() const{ return nIter; } 

template <typename T>  
double GradientConjugue<T>::getResVec1(int i) const{ return resVec1[i]; } 

template <typename T>  
double GradientConjugue<T>::getResVec(int i) const{ return resVec[i]; } 


template<typename T>
void GradientConjugue<T>::NormalisationResidu(){
	if(resolution){
		
		for(int i = 0; i<nIter; i++){
			resVec[i] = log10(resVec[i]/b.norm2());
		}
	}else{
		std::cout<<"Impossible de normaliser les residus car le probleme n a pas ete resolu"<<std::endl;
	}
}



#endif
