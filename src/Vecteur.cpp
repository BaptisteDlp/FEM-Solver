#include<iostream>
#include<vector>
#include<cassert>
#include<cstdlib>
#include<cmath>
#include"Vecteur.hpp"


using namespace std;



Vecteur::Vecteur(){
	m_taille = 1;
	m_coeff.push_back(0);
}



Vecteur::Vecteur(int taille){
	
	m_taille = taille;
	
	for(int i = 0; i<m_taille; i++){
		
		m_coeff.push_back(0);
		
	}
}



Vecteur::Vecteur(Vecteur const& v){
	
	m_taille = v.m_taille;
	
	for(int i = 0; i<v.m_taille; i++){
		m_coeff.push_back(v.m_coeff[i]);
	}
	
	
}



int Vecteur::getTaille() const{
	
	return m_taille;
	
	
}


void Vecteur::printVect() const{
	
	
	for(int i = 0; i<m_taille; i++){
			cout<<"|"<<m_coeff[i]<<"|"<<endl;
	}
	
	
}

void Vecteur::printFlux(std::ostream &flux) const
{
	flux<<endl;
	for(int i = 0; i<m_taille; i++){
			flux<<"|"<<m_coeff[i]<<"|"<<endl;
	}
	flux<<endl;
	
}


double& Vecteur::operator()(int i){
	
	assert(i>-1 && i<m_taille);
	return m_coeff[i];

}






double Vecteur::operator()(int i) const{
	
	assert(i>-1 && i<m_taille);
	return m_coeff[i];
	
	
}




Vecteur& Vecteur::operator+=(Vecteur const& v){
	
	assert(m_taille == v.m_taille);
	
	for(int i = 0; i<m_taille; i++){
		m_coeff[i] += v.m_coeff[i];
	}
	
	return *this;
}



Vecteur& Vecteur::operator-=(Vecteur const& v){
	
	assert(m_taille == v.m_taille);
	
	for(int i = 0; i<m_taille; i++){
		m_coeff[i] -= v.m_coeff[i];
	}
	
	return *this;
}




void Vecteur::ones(){
	
	for(int i = 0; i<m_taille; i++){
		m_coeff[i] = 1;
	}
	
}




double Vecteur::norm1() const{
	
	double result(0);
	
	for(int i = 0; i < m_taille; i++){
		result += abs(m_coeff[i]);
		
	}
	
	return result;
	
}

double Vecteur::norm2() const{
	
	double result(0);
	
	for(int i = 0; i < m_taille; i++){
		result += m_coeff[i]*m_coeff[i];
		
	}
	
	return sqrt(result);
	
}

double Vecteur::normInf() const{
	
	double result(abs(m_coeff[0]));
	
	for(int i = 0; i < m_taille-1; i++){
		if(result < abs(m_coeff[i+1])){
			result = abs(m_coeff[i+1]);
		}
		
	}
	
	return result;
	
}







Vecteur operator+(Vecteur const& v1, Vecteur const& v2){
	
	Vecteur a(v1);
	a+=v2;
	return a;
	
	
}




Vecteur operator-(Vecteur const& v1, Vecteur const& v2){
	
	Vecteur a(v1);
	a-=v2;
	return a;
	
	
}


Vecteur operator*(double c, Vecteur const& v){
	
	Vecteur a(v);
	
	for(int i = 0; i<a.getTaille(); i++){
		
			a(i) = c*a(i);
		
	}
	
	return a;
	
}

ostream& operator<<(std::ostream &flux, Vecteur const& v){
	
	v.printFlux(flux);
    return flux;
	
}







double operator*(Vecteur const& v1, Vecteur const& v2){
	
	double result(0);
	
	for(int i = 0; i<v1.getTaille(); i++){
		
			result += v1(i)*v2(i);
		
	}
	
	return result;
	
}




