#ifndef SPARSEMATRIXCOO2_HPP
#define SPARSEMATRIXCOO2_HPP

#include"Vecteur.hpp"

class SparseMatrixCoo2{
	
public:
  SparseMatrixCoo2();
  SparseMatrixCoo2(int const& taille);
  SparseMatrixCoo2(int nbR, int nbC, std::vector<int> row, std::vector<int> col, std::vector<double> val);
  SparseMatrixCoo2(SparseMatrixCoo2 const& m);
		
  void zeros();
  Vecteur MVProd(Vecteur const& v) const;
  Vecteur MVProdT(Vecteur const& v) const;
		
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
		
  void printRow() const;
  void printCol() const;
  void printVal() const;
		
  void printMat() const;
  void printFlux(std::ostream &flux) const;
		
		
		
		
		
		
		
		
		
	
private:
  int m_nbR;
  int m_nbC;
  int m_nnz;
  std::vector<int> m_row;
  std::vector<int> m_col;
  std::vector<double> m_val;
	
	
	
};



SparseMatrixCoo2 operator*(double c, SparseMatrixCoo2 const& m);
Vecteur operator*(SparseMatrixCoo2 const& m, Vecteur const& v);
std::ostream& operator<<(std::ostream &flux, SparseMatrixCoo2 const& m );





#endif
