#include<iostream>
#include<vector>
#include<cassert>
#include<stdlib.h>
#include"SparseMatrix.hpp"
#include"Vecteur.hpp"


using namespace std;




SparseMatrix::SparseMatrix(){
	
	m_nbR = 1;
	m_nbC = 1;
	m_nnz = 1;
	
	m_row.push_back(0);
	m_row.push_back(1);
	m_col.push_back(0);
	m_val.push_back(1.);
	
	
	
}
		
SparseMatrix::SparseMatrix(int const& taille){
	
	m_nbR = taille;
	m_nbC = taille;
	m_nnz = taille;
	
	//Changement ici
	for(int i = 0; i<taille; i++){
		m_row.push_back(i);
		m_col.push_back(i);
		m_val.push_back(1.);
		
		if(i == taille - 1){
			m_row.push_back(taille);
		}
	}
	
	
		
}





SparseMatrix::SparseMatrix(int nbR, int nbC, vector<int> row, vector<int> col, vector<double> val){
	
	
	m_nbR = nbR;
	m_nbC = nbC;
	
	
	for(unsigned int i = 0; i<row.size(); i++){
		
		m_row.push_back(row[i]);
		
	}
	
	for(unsigned int i = 0; i<col.size(); i++){
		
		m_col.push_back(col[i]);
		m_val.push_back(val[i]);
		
	}
	
	
	
	m_nnz = row[row.size()-1];
	
	
}
	
	
	
	

	
SparseMatrix::SparseMatrix(SparseMatrix const& m){
	
	m_nbR = m.getnbR();
	m_nbC = m.getnbC();
	m_nnz = m.getnnz();
	
	
	for(int i = 0; i<m.getnbR()+1; i++){
		
		m_row.push_back(m.getRow(i));
		
	}
	
	for(int i = 0; i<m.getnnz(); i++){
		
		m_col.push_back(m.getCol(i));
		m_val.push_back(m.getVal(i));
		
	}
	
	
	m_nnz = m.getRow(m.getRow(m.getnbR()));
	
	
}
	
	
	
	
	
	
void SparseMatrix::zeros(){
	
	for(unsigned int i = 0; i<m_val.size(); i++){
		m_val[i] = 0;
	}
	
}
	
	
	
	

	
Vecteur SparseMatrix::MVProd(Vecteur const& v){
	
	assert(m_nbR == m_nbC);
	assert(m_nbC == v.getTaille());
	
	Vecteur a(v.getTaille());
	
	int k = 0;
	
	for(int i = 0; i<v.getTaille(); i++){
		
		for(int j = m_row[i]; j<m_row[i+1]; j++){
			a(i) += m_val[k]*v(m_col[j]);
			k++;
		}
		
		
	}
	
	return a;
	
	
	
}
	
	
	
	
	


int SparseMatrix::getnbR() const{ return m_nbR;}

int SparseMatrix::getnbC() const{ return m_nbC;}

int SparseMatrix::getnnz() const{ return m_nnz;}

int SparseMatrix::getRow(int i) const{ return m_row[i];}

int SparseMatrix::getCol(int i) const{ return m_col[i];}

double SparseMatrix::getVal(int i) const{ return m_val[i];}






void SparseMatrix::printRow() const{
	
	cout<<"Row:";
	for( unsigned int i = 0; i<m_row.size(); i++){
		cout<<"|"<<m_row[i];
	}
	cout<<"|"<<endl;
	
}

void SparseMatrix::printCol() const{
	
	cout<<"Col:";
	for(unsigned int i = 0; i<m_col.size(); i++){
		cout<<"|"<<m_col[i];
	}
	cout<<"|"<<endl;
	
}


void SparseMatrix::printVal() const{
	
	cout<<"Val:";
	for(unsigned int i = 0; i<m_val.size(); i++){
		cout<<"|"<<m_val[i];
	}
	cout<<"|"<<endl;
	
}






		
void SparseMatrix::printMat() const{
	
	
	for(int i = 0; i<m_nbR; i++){
		cout<<"|";
		for(int j = 0; j<m_nbC; j++){
			
			cout<<" "<<(*this)(i,j)<<" ";
			
		}
		cout<<"|"<<endl;
	}
	
	
}





		
		
void SparseMatrix::printFlux(ostream &flux) const
{
	flux<<endl;
    for(int i = 0; i<m_nbR; i++){
		flux<<"| ";
		for(int j = 0; j<m_nbC; j++){
			
			flux<<" "<<(*this)(i,j)<<" ";
			
		}
		flux<<" |"<<endl;
	}
	flux<<endl;
}





		
double SparseMatrix::operator()(int i, int j) const{
	
	assert(i>-1 && i<m_nbR && j>-1 && j<m_nbC);
	
	int k(0), l(0), m(0);
	bool ok = false;
	
	k = m_row[i]; 
	l = m_row[i+1]; //???
	
	for(m = k; m<l; m++){
		
		if(m_col[m] == j){
			ok = true;
			break;
		}
	}
	
	if(ok){
		return m_val[m];
	}else{
		return 0.;
	}
	
	
}
		
void SparseMatrix::setVal(int i, int j, double val){
	
	assert(i>-1 && i<m_nbR && j>-1 && j<m_nbC);
	
	int k(0), l(0), m(0);
	bool ok = false;
	
	k = m_row[i]; 
	l = m_row[i+1]; //???
	
	for(m = k; m<l; m++){
		
		if(m_col[m] == j){
			ok = true;
			break;
		}
	}
	
	if(ok){
		m_val[m] += val;
	}else{
		exit(0);
	}
	
}




		
		
Vecteur operator*(SparseMatrix const& m, Vecteur const& v){
	
	assert(m.getnbC() == v.getTaille());
	
	Vecteur a(v.getTaille());
	
	int k = 0;
	
	for(int i = 0; i<v.getTaille(); i++){
		
		for(int j = m.getRow(i); j<m.getRow(i+1); j++){
			a(i) += m.getVal(k)*v(m.getCol(j));
			k++;
		}
		
		
	}
	
	return a;
	

	
}
		
		
		
		




		
ostream& operator<<( ostream &flux, SparseMatrix const& m ){
    m.printFlux(flux);
    return flux;
}		
		
		





		
		
Vecteur inv_triang_sup(SparseMatrix const& m, Vecteur const& y){
	
	assert(m.getnbR() == m.getnbC());
	assert(m.getnbC() == y.getTaille());
	
	Vecteur z(y);
	
	
	for(int i = z.getTaille()-1; i>-1; i--){
		
		z(i) = y(i)/m(i,i);
		
		for(int j = z.getTaille()-1; j>i; j--){
			z(i) -= (m(i,j)*z(j))/m(i,i);
		}
		
	}
	
	return z;
}
		
		
		
		
		
		



		
		
SparseMatrix matriceLaplacienc(int taille){
	
	assert(taille > 1);
	
	vector<int> row;
	vector<int> col;
	vector<double> val;
	
	int nnz = 3*taille -2;
	
	row.push_back(0);
	row.push_back(2);
	
	col.push_back(0);
	col.push_back(1);
	
	val.push_back(2);
	val.push_back(-1);
	
	int i = nnz - 2;
	int column = 2;
	int r = 5;
	
	for(int i=2; i<taille+1; i++){
		row.push_back(r);
		if(i != taille -1){
			r += 3; 
		}else{
			r += 2;
		}
	}
	
	while(i>2){
		col.push_back(column - 2);
		col.push_back(column - 1);
		col.push_back(column);
		val.push_back(-1);
		val.push_back(2);
		val.push_back(-1);
		column += 1;
		i -= 3;
		
	}
	
	
	
	col.push_back(taille-2);
	col.push_back(taille-1);
	val.push_back(-1);
	val.push_back(2);
	
	return SparseMatrix(taille,taille,row,col,val);
	
	
}
		
		
		
		
		
		

		
		
		
		
SparseMatrix matalea(double p, long long unsigned int taille){
	
	
	assert(taille > 1 && p>0 && p<=1);
	
	vector<int> row;
	vector<int> col;
	vector<double> val;
	
	
	long long unsigned int nnz = p*taille*taille;
	
	
	long long unsigned int nbrow = nnz / taille;
	
	
	long long unsigned int k = nbrow;
	long long unsigned int place = taille / nbrow;
	long long unsigned int l(0);
	
	
	
	long long unsigned int i = 0;
	
	row.push_back(0);
	
	
	
	for(long long unsigned int i=1; i<taille+1; i++){
		
		row.push_back(k);
		k += nbrow;
		
	}
	
	
	while(i<nnz){
		
		l = 0;
		
		for(long long unsigned int j = 0; j<nbrow; j++){
			
			col.push_back(l);
			val.push_back( rand() % 100 +1);
			l += place;
			i++;
		}

    }
	
	return SparseMatrix(taille,taille,row,col,val);
	
}







		
	

SparseMatrix matHnc(int taille){
	
	
	assert(taille > 1);
	
	vector<int> row;
	vector<int> col;
	vector<double> val;
	int nbcoeff(0);
	
	row.push_back(0);
	
	for(int i=0; i<taille; i++){
		nbcoeff = 0;
		
		if(i-2>=0){
			col.push_back(i-2);
			val.push_back(1);
			nbcoeff++;
		}
		
		col.push_back(i);
		val.push_back(i+1);
		nbcoeff++;
		
		if(i+2<taille){
			col.push_back(i+2);
			val.push_back(1);
			nbcoeff++;
		}
		
		row.push_back(row[i]+nbcoeff);
		
	}
	
	return SparseMatrix(taille,taille,row,col,val);
	

	
}
		
	

	
	
	
	
	
	
	

SparseMatrix PrecondJacobiC(SparseMatrix const& m){

	assert(m.getnbR() == m.getnbC());
	assert(m.getnbR() > 1);
	
	vector<int> row;
	vector<int> col;
	vector<double> val;
	
	for(int i = 0; i<m.getnbR(); i++){
		row.push_back(i);
		col.push_back(i);
		val.push_back(1/m(i,i));
		
		if(i == m.getnbR() - 1){
			row.push_back(m.getnbR());
		}
		
	}

	return SparseMatrix(m.getnbR(),m.getnbC(),row,col,val);

}	
		
		

	