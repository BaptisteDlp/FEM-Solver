#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include"DenseMatrix.hpp"
#include"SparseMat.hpp"
#include"Vecteur.hpp"
#include"mesh.hpp"


double areaTri(double node1X,double node1Y,double node2X,double node2Y,double node3X,double node3Y);

DenseMatrix MassElem(double node1X,double node1Y,double node2X,double node2Y,double node3X,double node3Y);
DenseMatrix RigElem(double node1X,double node1Y,double node2X,double node2Y,double node3X,double node3Y);

MatCSR Stiffness_Mass_matrix(Mesh const& Th, bool stiff = true);







Vecteur setSecondMember(DenseMatrix const& tabNodes);
Vecteur AssemblyF(DenseMatrix const& M, Vecteur const& f);
double quadratureTri(int i,double node1X, double node1Y, double node2X, double node2Y, double node3X, double node3Y);
Vecteur assemblyFQuad(DenseMatrix tabNodes, DenseMatrix tabElts);
DenseMatrix quadratureTriG(double node1X, double node1Y, double node2X, double node2Y, double node3X, double node3Y);
DenseMatrix AssemblyMG(DenseMatrix const& tabNodes,DenseMatrix const& tabElts);


#endif
