#ifndef MESH_HPP
#define MESH_HPP


#include"DenseMatrix.hpp"
#include"Vecteur.hpp"
#include<utility>

DenseMatrix loadNodes(std::string fileName, bool offset = false);
DenseMatrix loadElements(std::string fileName, bool offset = false);


Vecteur setX(DenseMatrix const& tabNodes);
Vecteur setY(DenseMatrix const& tabNodes);

void rectMesh(std::string fileName, int nbNodesH, int nbNodesV, double x0, double xN, double y0, double yN);


std::vector<std::pair<int,std::pair<double,double>>> meshExtraction(std::string const&  meshNameFile ,std::string const& regionNameFile, std::string const& meshLocNameFile, int regionLabelIn, int regionLabelOut, int labelIn, int labelOut, int labelCommon = -1);

std::vector<std::pair<int,int>> buildGammaij(DenseMatrix const& tabv1, DenseMatrix const& tabv2,int labelIn);

std::vector<std::pair<int,int>> buildGamma(std::vector<std::pair<int,std::pair<double,double>>> const& g, DenseMatrix const& tabv1, DenseMatrix const& tabv2 ,int labelIn);


void exportMeshData(std::string meshName, DenseMatrix const& tabNodes, DenseMatrix const& tabElts);
void exportMeshGnuplot(std::string scriptName, std::string meshName);

void exportSolData(std::string SolName, DenseMatrix const& tabNodes, DenseMatrix const& tabElts, Vecteur const& solution);
void exportSolGnuplot(std::string scriptName, std::string SolName, std::string title);



std::vector<std::vector<double>> sol(DenseMatrix const& tabNodes, Vecteur const& solution);
void exportSolDataPm3d(std::string SolName, DenseMatrix const& tabNodes, Vecteur const& solution);
void exportSolGnuplotPm3d(std::string scriptName, std::string SolName, std::string title);







#endif
