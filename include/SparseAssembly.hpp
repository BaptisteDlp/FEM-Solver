#ifndef SPARSEASSEMBLY_HPP
#define SPARSEASSEMBLY__HPP

#include<vector>
#include"SparseMatrixCoo.hpp"
#include"SparseMatrixCoo2.hpp"
#include"SparseMatrix.hpp"
#include"DenseMatrix.hpp"
#include"Vecteur.hpp"



void mergeRow(std::vector<int> &row, std::vector<int> &col, std::vector<double> &val, int l, int m, int r);
void mergeSortRow(std::vector<int> &row, std::vector<int> &col, std::vector<double> &val, int l, int r);


SparseMatrixCoo2 sparseMatrixCoo2NeumannH(DenseMatrix const& tabNodes, DenseMatrix const& tabElts);


SparseMatrixCoo sparseMatrixCooDirichletH(DenseMatrix const& tabNodes, DenseMatrix const& tabElts, std::vector<int> const& boundary);
SparseMatrixCoo2 sparseMatrixCoo2Laplace(DenseMatrix const& tabNodes, DenseMatrix const& tabElts, std::vector<int> const& boundary);

SparseMatrix sparseMatrixDirichletH(DenseMatrix const& tabNodes, DenseMatrix const& tabElts, std::vector<Vecteur> const& tabEdges, std::vector<int> const& boundary);




#endif