#ifndef SPARSEMATRIXCOO_HPP
#define SPARSEMATRIXCOO_HPP

#include"Vecteur.hpp"

class SparseMatrixCoo{
	
public:
  SparseMatrixCoo();
  SparseMatrixCoo(int const& taille);
  SparseMatrixCoo(int nbR, int nbC, std::vector<int> row, std::vector<int> col, std::vector<double> val);
  SparseMatrixCoo(SparseMatrixCoo const& m);
		
  void zeros();
  Vecteur MVProd(Vecteur const& v) const;
		
  int getnbR() const;
  int getnbC() const;
  int getnnz() const;
  int getRow(int i) const;
  int getCol(int i) const;
  double getVal(int i) const;
  std::vector<int> get_row() const;
  std::vector<int> get_col() const;
  std::vector<double> get_val() const;
		
		
  void Insert(int row, int col, double val);
  void add(int i, double val){
    m_val[i] += val;
  }
		
  void printRow() const;
  void printCol() const;
  void printVal() const;
		
  void printMat() const;
  void printFlux(std::ostream &flux) const;
		
		
  double operator()(int i, int j) const;
		
		
		
		
		
		
	
private:
  int m_nbR;
  int m_nbC;
  int m_nnz;
  std::vector<int> m_row;
  std::vector<int> m_col;
  std::vector<double> m_val;
	
	
	
};



SparseMatrixCoo operator*(double c, SparseMatrixCoo const& m);
Vecteur operator*(SparseMatrixCoo const& m, Vecteur const& v);
std::ostream& operator<<(std::ostream &flux, SparseMatrixCoo const& m );

Vecteur inv_triang_sup(SparseMatrixCoo const& m, Vecteur const& y);
SparseMatrixCoo matriceLaplaciencoo(int taille);






#endif
