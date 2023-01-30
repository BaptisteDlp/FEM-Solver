#ifndef SURFACICASSEMBLY_HPP
#define SURFACICASSEMBLY_HPP

#include<vector>
#include"DenseMatrix.hpp"
#include"Vecteur.hpp"


double lengthEdge(double node1X, double node1Y, double node2X , double node2Y);
Vecteur boundElemLin(double node1X, double node1Y, double node2X , double node2Y);
Vecteur assemblyL(DenseMatrix const& tabNodes, std::vector<Vecteur> const& tabEdges);
DenseMatrix boundElemBil(double node1X, double node1Y, double node2X , double node2Y);
DenseMatrix assemblageBil(DenseMatrix const& tabNodes, std::vector<Vecteur> const& tabBound);


#endif