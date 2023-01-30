#ifndef VECTEUR_HPP
#define VECTEUR_HPP


#include<vector>




class Vecteur{
	
	public:
		Vecteur();
		Vecteur(int taille);
		Vecteur(Vecteur const& m);
		
		int getTaille() const;
		void printVect() const;
		void printFlux(std::ostream &flux) const;
		
		double& operator()(int i);
		double operator()(int i) const;
		
		Vecteur& operator+=(Vecteur const& v);
		Vecteur& operator-=(Vecteur const& v);
		
		void ones(); 
		double norm1() const;
		double norm2() const;
		double normInf() const;
		
	
	private:
		int m_taille;
		std::vector<double> m_coeff;
	
	
	
};

Vecteur operator+(Vecteur const& v1, Vecteur const& v2);
Vecteur operator-(Vecteur const& v1, Vecteur const& v2);
Vecteur operator*(double c, Vecteur const& v);
double operator*(Vecteur const& v1, Vecteur const& v2);
std::ostream& operator<<(std::ostream &flux, Vecteur const& m );



#endif