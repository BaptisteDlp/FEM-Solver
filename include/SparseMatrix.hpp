#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include"Vecteur.hpp"

class SparseMatrix{
	
	public:
		SparseMatrix();
		SparseMatrix(int const& taille);
		SparseMatrix(int nbR, int nbC, std::vector<int> row, std::vector<int> col, std::vector<double> val);
		SparseMatrix(SparseMatrix const& m);
		
		void zeros();
		Vecteur MVProd(Vecteur const& v);
		
		int getnbR() const;
		int getnbC() const;
		int getnnz() const;
		int getRow(int i) const;
		int getCol(int i) const;
		double getVal(int i) const;
		
		
		void printRow() const;
		void printCol() const;
		void printVal() const;
		
		void printMat() const;
		void printFlux(std::ostream &flux) const;
		
		
		double operator()(int i, int j) const;
		void setVal(int i, int j, double val);
		
		
		
		
		
		
	
	private:
		int m_nbR;
		int m_nbC;
		int m_nnz;
		std::vector<int> m_row;
		std::vector<int> m_col;
		std::vector<double> m_val;
	
	
	
};


Vecteur operator*(SparseMatrix const& m, Vecteur const& v);
SparseMatrix operator*(double c, SparseMatrix const& m);
std::ostream& operator<<(std::ostream &flux, SparseMatrix const& m );

Vecteur inv_triang_sup(SparseMatrix const& m, Vecteur const& y);
SparseMatrix matriceLaplacienc(int taille);
SparseMatrix matalea(double p, long long unsigned int taille);
SparseMatrix matHnc(int taille);
SparseMatrix PrecondJacobiC(SparseMatrix const& m);





#endif