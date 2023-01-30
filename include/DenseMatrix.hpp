#ifndef DENSEMATRIX_HPP
#define DENSEMATRIX_HPP



#include<iostream>
#include"Vecteur.hpp"


class DenseMatrix{
	
public:
  DenseMatrix();
  DenseMatrix(int nbR, int nbC);
  DenseMatrix(DenseMatrix const& m);
		
  void Load(char* const filename);
		
  int getnbR() const;
  int getnbC() const;
  double getCoeff(int i) const;
  bool getFactoLU() const;
  bool getFactoCho() const;
		
		
  void printMat() const;
  void print_nnz() const;
  void printFlux(std::ostream &flux) const;
		
  void decomp_LU();
  void decomp_Cholesky();
		
  double& operator()(int i, int j);
  double operator()(int i, int j) const;
		
  DenseMatrix& operator+=(DenseMatrix const& m);
  DenseMatrix& operator-=(DenseMatrix const& m);
		
  Vecteur MVProdT(Vecteur const& v);
		
	
private:
  int m_nbR;
  int m_nbC;
  std::vector<double> m_coeff;
  bool m_factoLU;
  bool m_factoCho;
		
  friend double FrobNorm(DenseMatrix const& M);
	
	
	
};

DenseMatrix operator+(DenseMatrix const& m1, DenseMatrix const& m2);
DenseMatrix operator-(DenseMatrix const& m1, DenseMatrix const& m2);
DenseMatrix operator*(double c, DenseMatrix const& m);
DenseMatrix operator*(DenseMatrix const& m1, DenseMatrix const& m2);
std::ostream& operator<<(std::ostream &flux, DenseMatrix const& m );

Vecteur inv_triang_inf (DenseMatrix const& m, Vecteur const& y);
Vecteur inv_triang_sup(DenseMatrix const& m, Vecteur const& y);

Vecteur inv_LU(DenseMatrix &a, Vecteur const& b);
Vecteur inv_Cho(DenseMatrix &a, Vecteur const& b);

void inv_LU(DenseMatrix &a, double *b, double *z);
void inv_Cho(DenseMatrix &a, double *b, double *z);

Vecteur solveLU(DenseMatrix &a, Vecteur const& b);
Vecteur solveCholesky(DenseMatrix &a, Vecteur const& b);

void solveLU(DenseMatrix &a, double* b, double *z);
void solveCholesky(DenseMatrix &a, double *b, double *z);

Vecteur operator*(DenseMatrix const& m, Vecteur const& v);
double* operator*(DenseMatrix const& m, double * v);

DenseMatrix matLaplacienD(int taille);
DenseMatrix matHn(int taille);
DenseMatrix PrecondJacobi(DenseMatrix const& A);

DenseMatrix B_Matrix(int size);




#endif
