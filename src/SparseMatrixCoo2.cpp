#include<iostream>
#include<vector>
#include<cassert>
#include<stdlib.h>
#include"SparseMatrixCoo2.hpp"
#include"Vecteur.hpp"


using namespace std;




SparseMatrixCoo2::SparseMatrixCoo2(){
	
	m_nbR = 1;
	m_nbC = 1;
	m_nnz = 1;
	
	m_row.push_back(0);
	m_col.push_back(0);
	m_val.push_back(1.);
	
	
	
}
	


	
SparseMatrixCoo2::SparseMatrixCoo2(int const& taille){
	
	m_nbR = taille;
	m_nbC = taille;
	m_nnz = taille;
	

	for(int i = 0; i<taille; i++){
		m_row.push_back(i);
		m_col.push_back(i);
		m_val.push_back(1.);
	}
	
	
		
}





SparseMatrixCoo2::SparseMatrixCoo2(int nbR, int nbC, vector<int> row, vector<int> col, vector<double> val){
	
	
	m_nbR = nbR;
	m_nbC = nbC;
	
	
	
	
	for(unsigned int i = 0; i<col.size(); i++){
		
		m_row.push_back(row[i]);
		m_col.push_back(col[i]);
		m_val.push_back(val[i]);
		
	}
	
	
	
	m_nnz = val.size();
	
	
}
	
	
	
	

	
SparseMatrixCoo2::SparseMatrixCoo2(SparseMatrixCoo2 const& m){
	
	m_nbR = m.getnbR(); 
	m_nbC = m.getnbC(); 
	m_nnz = m.getnnz(); 
	
	for(int i = 0; i<m.getnnz(); i++){
		m_row.push_back(m.getRow(i)); 
		m_col.push_back(m.getCol(i));
		m_val.push_back(m.getVal(i));
	}
	
	
	
	
}
	
	
	
	
	
	
void SparseMatrixCoo2::zeros(){
	
	vector<int>::iterator it2 = m_col.begin();
	vector<double>::iterator it3 = m_val.begin();
	
	for(vector<int>::iterator it1 = m_row.begin(); it1 != m_row.end(); ++it1){
		m_row.erase(it1);
		m_col.erase(it2);
		m_val.erase(it3);
		it2++;
		it3++;
	}
	
}
	
	
	
	

	
Vecteur SparseMatrixCoo2::MVProd(Vecteur const& v) const{
	
	assert(m_nbR == m_nbC);
	assert(m_nbC == v.getTaille());
	
	Vecteur a(v.getTaille());
	
	
	
	for(unsigned l = 0; l<m_val.size(); l++){
		
		a(m_row[l]) += m_val[l]*v(m_col[l]);
		
	}
	
	
	return a;
	
}
	
	
	
	
	


int SparseMatrixCoo2::getnbR() const{ return m_nbR;}

int SparseMatrixCoo2::getnbC() const{ return m_nbC;}

int SparseMatrixCoo2::getnnz() const{ return m_nnz;}

int SparseMatrixCoo2::getRow(int i) const{ return m_row[i];}

int SparseMatrixCoo2::getCol(int i) const{ return m_col[i];}

double SparseMatrixCoo2::getVal(int i) const{ return m_val[i];}

std::vector<int> SparseMatrixCoo2::get_row() const{ return m_row;}
std::vector<int> SparseMatrixCoo2::get_col() const{ return m_col;}
std::vector<double> SparseMatrixCoo2::get_val() const{ return m_val;}





void SparseMatrixCoo2::Insert(int row, int col, double val){
	
	assert(row>-1 && row<m_nbR && col>-1 && col<m_nbC && val != 0);
	

	
	m_row.push_back(row);
	m_col.push_back(col);
	m_val.push_back(val);
	m_nnz++;
	
	
	
}








void SparseMatrixCoo2::printRow() const{
	
	cout<<"Row:";
	for( unsigned int i = 0; i<m_row.size(); i++){
		cout<<"|"<<m_row[i];
	}
	cout<<"|"<<endl;
	
}

void SparseMatrixCoo2::printCol() const{
	
	cout<<"Col:";
	for(unsigned int i = 0; i<m_col.size(); i++){
		cout<<"|"<<m_col[i];
	}
	cout<<"|"<<endl;
	
}


void SparseMatrixCoo2::printVal() const{
	
	cout<<"Val:";
	for(unsigned int i = 0; i<m_val.size(); i++){
		cout<<"|"<<m_val[i];
	}
	cout<<"|"<<endl;
	
}












Vecteur operator*(SparseMatrixCoo2 const& m, Vecteur const& v){
	
	return m.MVProd(v);
	
}

		
		
		
Vecteur SparseMatrixCoo2::MVProdT(Vecteur const& v) const{
	
	assert(m_nbR == m_nbC);
	assert(m_nbC == v.getTaille());
	
	Vecteur a(v.getTaille());
	
	
	
	for(unsigned l = 0; l<m_val.size(); l++){
		
		a(m_col[l]) += m_val[l]*v(m_row[l]);
		
	}
	
	
	return a;
	
}

		
		
		
		




	
		
		



		
		
		

		
		






	
		
	

	
	
	
	
	
	
	
		

	