#ifndef EDGE_HPP
#define EDGE_HPP

#include<vector>
#include"DenseMatrix.hpp"
#include"Vecteur.hpp"


std::vector<Vecteur> findEdge(DenseMatrix const& tabNodes, DenseMatrix const& tabElts);
std::vector<Vecteur> findBoundary(DenseMatrix const& tabNodes, DenseMatrix const& tabElts);
std::vector<int> findNodesBoundary(DenseMatrix const& tabNodes, std::vector<Vecteur> const& tabBoundary);
std::vector<int> idBord(DenseMatrix const& tabNodes, std::vector<int> const& tabNodesBound);

void printEdges(std::vector<Vecteur> const& edge);
void printBound(std::vector<Vecteur> const& bound);
void printNodesBound(std::vector<int> const& nodesBound);
void printIdBord(std::vector<int> const& idBord);




#endif