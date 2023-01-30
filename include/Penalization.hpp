#ifndef PENALIZATION_HPP
#define PENALIZATION_HPP

#include<vector>
#include"DenseMatrix.hpp"
#include"Vecteur.hpp"


void pseudoElemA(DenseMatrix &A, std::vector<int> const& boundary);
void penalizationF(Vecteur &F, std::vector<int> const& boundary);
void penalizationA(DenseMatrix &A, std::vector<int> const& boundary);
void penalizationFNonHomo(Vecteur &F, std::vector<int> const& boundary,DenseMatrix const& tabNodes);



DenseMatrix periodicConditionA(DenseMatrix const& A, std::vector<int> const& tabIdBord);
void periodicConditionA2(DenseMatrix &A, std::vector<int> const& tabIdBord);
void periodicConditionF(Vecteur &F, std::vector<int> const& tabIdBord);
void periodicConditionSol(Vecteur &Sol, std::vector<int> const &tabIdBord);






#endif