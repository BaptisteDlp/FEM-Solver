#include<iostream>
#include<cmath>
#include"SparseAssembly.hpp"
#include"SurfacicAssembly.hpp"
#include"Assembly.hpp"


using namespace std;



void mergeRow(vector<int> &row, vector<int> &col, vector<double> &val, int l, int m, int r){
  unsigned int i, j, k;
  int n1 = m - l + 1;
  int n2 =  r - m;
 
  /* create temp arrays */
  vector<int> L1(n1,0), R1(n2,0), L2(n1,0), R2(n2,0);
  vector<double> L3(n1,0), R3(n2,0);
	
  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++){
    L1[i] = row[l + i];
    L2[i] = col[l + i];
    L3[i] = val[l + i];
  }
  for (j = 0; j < n2; j++){
    R1[j] = row[m + 1+ j];
    R2[j] = col[m + 1+ j];
    R3[j] = val[m + 1+ j];
  }
 
  /* Merge the temp arrays back into arr[l..r]*/
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while (i < n1 && j < n2)
    {
      if (L1[i] <= R1[j])
        {
	  row[k] = L1[i];
	  col[k] = L2[i];
	  val[k] = L3[i];
	  i++;
        }
      else
        {
	  row[k] = R1[j];
	  col[k] = R2[j];
	  val[k] = R3[j];
	  j++;
        }
      k++;
    }
 
  /* Copy the remaining elements of L[], if there
     are any */
  while (i < n1)
    {
      row[k] = L1[i];
      col[k] = L2[i];
      val[k] = L3[i];
      i++;
      k++;
    }
 
  /* Copy the remaining elements of R[], if there
     are any */
  while (j < n2)
    {
      row[k] = R1[j];
      col[k] = R2[j];
      val[k] = R3[j];
      j++;
      k++;
    }
	
}
 
 
 
 
void mergeSortRow(vector<int> &row, vector<int> &col, vector<double> &val, int l, int r){
	
  if (l < r){
    int m = l+(r-l)/2; 
 
    mergeSortRow(row,col,val, l, m);
    mergeSortRow(row,col,val, m+1, r);
 
    mergeRow(row,col,val, l, m, r);
  }
	
}









SparseMatrixCoo2 sparseMatrixCoo2NeumannH(DenseMatrix const& tabNodes, DenseMatrix const& tabElts){
	
	
  vector<int> row;
  vector<int> col;
  vector<double> val;
	
  SparseMatrixCoo2 A(tabNodes.getnbR(),tabNodes.getnbR(),row,col,val);

   
  for(int i = 0; i<tabElts.getnbR(); i++){
    //cout<<i<<endl;
    double node1X = tabNodes(tabElts(i,0),0);
    double node1Y = tabNodes(tabElts(i,0),1);
    double node2X = tabNodes(tabElts(i,1),0);
    double node2Y = tabNodes(tabElts(i,1),1);
    double node3X = tabNodes(tabElts(i,2),0);
    double node3Y = tabNodes(tabElts(i,2),1);
        
    DenseMatrix Mt = MassElem(node1X,node1Y,node2X,node2Y,node3X,node3Y);
    DenseMatrix Kt = RigElem(node1X,node1Y,node2X,node2Y,node3X,node3Y);

    for(int j = 0; j<3; j++){
      for(int k = 0; k<3; k++){ 
	double p = Mt(j,k) + Kt(j,k);
	if(p != 0){
	  A.Insert(tabElts(i,j),tabElts(i,k),p);
	}			
      }
    }
  }
	
	
	
  return A;
}












SparseMatrixCoo sparseMatrixCooDirichletH(DenseMatrix const& tabNodes, DenseMatrix const& tabElts, vector<int> const& boundary){
	
	
  vector<int> row;
  vector<int> col;
  vector<double> val;
	
  SparseMatrixCoo A(tabNodes.getnbR(),tabNodes.getnbR(),row,col,val);

   
  for(int i = 0; i<tabElts.getnbR(); i++){
    //cout<<i<<endl;
    double node1X = tabNodes(tabElts(i,0),0);
    double node1Y = tabNodes(tabElts(i,0),1);
    double node2X = tabNodes(tabElts(i,1),0);
    double node2Y = tabNodes(tabElts(i,1),1);
    double node3X = tabNodes(tabElts(i,2),0);
    double node3Y = tabNodes(tabElts(i,2),1);
        
    DenseMatrix Mt = MassElem(node1X,node1Y,node2X,node2Y,node3X,node3Y);
    DenseMatrix Kt = RigElem(node1X,node1Y,node2X,node2Y,node3X,node3Y);

    for(int j = 0; j<3; j++){
      for(int k = 0; k<3; k++){ 
	double p = Mt(j,k) + Kt(j,k);
	if(p != 0){
	  A.Insert(tabElts(i,j),tabElts(i,k),p);
	}		
      }
    }
  }
	
  //Penalisation
  for(unsigned int j = 0; j<boundary.size(); j++){
    double p = pow(10,30)-A(boundary[j],boundary[j]);
    A.Insert(boundary[j],boundary[j],p);	
  }		
  
	
  return A;
}




SparseMatrixCoo2 sparseMatrixCoo2Laplace(DenseMatrix const& tabNodes, DenseMatrix const& tabElts, vector<int> const& boundary){
	
	
  vector<int> row;
  vector<int> col;
  vector<double> val;
	
  SparseMatrixCoo2 A(tabNodes.getnbR(),tabNodes.getnbR(),row,col,val);

   
  for(int i = 0; i<tabElts.getnbR(); i++){
    double node1X = tabNodes(tabElts(i,0),0);
    double node1Y = tabNodes(tabElts(i,0),1);
    double node2X = tabNodes(tabElts(i,1),0);
    double node2Y = tabNodes(tabElts(i,1),1);
    double node3X = tabNodes(tabElts(i,2),0);
    double node3Y = tabNodes(tabElts(i,2),1);
        
    DenseMatrix Kt = RigElem(node1X,node1Y,node2X,node2Y,node3X,node3Y);

    for(int j = 0; j<3; j++){
      for(int k = 0; k<3; k++){ 
	double p =  Kt(j,k);
	if(p != 0){
	  A.Insert(tabElts(i,j),tabElts(i,k),p);
	}			
      }
    }
  }
	
  //Penalisation
  for(unsigned int j = 0; j<boundary.size(); j++){
    double p = pow(10,30);
    A.Insert(boundary[j],boundary[j],p);	
  }		
  
	
  return A;
}
		

		
		
		
		
		
		
		
		

SparseMatrix sparseMatrixDirichletH(DenseMatrix const& tabNodes, DenseMatrix const& tabElts, vector<Vecteur> const& tabEdges, vector<int> const& boundary){
	
	
  vector<int> row;
  vector<int> col;
  vector<double> val;
	
  for(unsigned int i = 0; i<tabEdges.size(); i++){
		
    row.push_back(tabEdges[i](0));
    col.push_back(tabEdges[i](1));
    val.push_back(0.);

    row.push_back(tabEdges[i](1));
    col.push_back(tabEdges[i](0));
    val.push_back(0);
  }
	
  for(int i = 0; i<tabNodes.getnbR(); i++){
		
    row.push_back(i);
    col.push_back(i);
    val.push_back(0);
  }
	
  mergeSortRow(row,col,val,0,row.size()-1);
	
  unsigned i(0);
	
  while(i<row.size()){
    int first = row[i]; int begin(i);
    while(i<row.size() && row[i] == first){
      i++;
    }
    mergeSortRow(col,row,val,begin,i-1);
  }
	
  i = 1;
	
	
  vector<int>::iterator it = row.begin();
  vector<int> newRow(tabNodes.getnbR()+1,0);
  int first = 0; int l = 0; int nb(0);

  while(it != row.end()){
    l = *it;
    if(l == first){
      nb++; it++; 
    }else{
      newRow[i] = nb; first++; i++; 
    }
		
  }
  newRow[i] = nb;
	
	
  SparseMatrix A(tabNodes.getnbR(),tabNodes.getnbR(),newRow,col,val);
	
	
   
  for(int i = 0; i<tabElts.getnbR(); i++){
    double node1X = tabNodes(tabElts(i,0),0);
    double node1Y = tabNodes(tabElts(i,0),1);
    double node2X = tabNodes(tabElts(i,1),0);
    double node2Y = tabNodes(tabElts(i,1),1);
    double node3X = tabNodes(tabElts(i,2),0);
    double node3Y = tabNodes(tabElts(i,2),1);
        
    DenseMatrix Mt = MassElem(node1X,node1Y,node2X,node2Y,node3X,node3Y);
    DenseMatrix Kt = RigElem(node1X,node1Y,node2X,node2Y,node3X,node3Y);
		
    for(int j = 0; j<3; j++){
      for(int k = 0; k<3; k++){ 
	double p = Mt(j,k) + Kt(j,k);
	A.setVal(tabElts(i,j),tabElts(i,k),p);
      }
    }
  }
	
	
  cout<<A<<endl;
  A.printRow(); A.printCol(); A.printVal();
	
  for(unsigned int j = 0; j<boundary.size(); j++){
    double p = pow(10,30)-A(boundary[j],boundary[j]);
    A.setVal(boundary[j],boundary[j],p);	
  }	
		
	
	
  return A;
}







