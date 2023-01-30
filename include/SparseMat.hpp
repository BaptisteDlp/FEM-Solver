#include "Vecteur.hpp"
#include <map>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

extern int debug ;



struct  MatVirt;

struct MatVirt {
public:
    int n,m;
  virtual  Vecteur  addmatmul(Vecteur const& x, Vecteur & Ax) const =0;
    virtual ~MatVirt() {}

    Vecteur matmul(Vecteur const& x, Vecteur & Ax) const
    {
      for( int i = 0; i<n; i++){
	Ax(i) = 0;
      }
      return addmatmul(x,Ax);
    }
    MatVirt(int nn,int mm=-1) : n(nn),m(mm<0 ?nn:mm) {}
};



typedef std::map<std::pair<int,int>,double> MatrixMap;



struct MatCSR : public MatVirt
{
    int nnz;
    std::vector<int> ip,j;// ip pointeur de debut de ligne , j n* col
    std::vector<double> a; // coef
    MatCSR(int nn,int mm,const MatrixMap &A) :MatVirt(nn,mm),nnz(A.size()),ip(n+1),j(nnz),a(nnz)
    { // c++11  -std=c++11
        int k=0;
        std::fill(ip.begin(), ip.end(), -1);// unset
        ip[0] =0;
        for (auto l=A.begin(); l != A.end(); ++l )
        {
            int i = l->first.first;
            j[k] = l->first.second;
            a[k] = l->second;
            ip[i+1] = ++k;
        }
        // nettoyage
        for (int i=0; i<n; ++i)
            if(ip[i+1] ==-1) // not unset ...
                ip[i+1] = ip[i];
        
        assert(ip[n]==nnz);
        std::cout << " MatCSR " << n << " " << m << " " << nnz << std::endl;
    }

    MatCSR(MatCSR const& A): MatVirt(A.n,A.m),nnz(A.nnz),ip(A.n+1),j(A.nnz),a(A.nnz){
      for(unsigned int i = 0; i<ip.size(); i++){
	ip[i] = A.ip[i];
      }

      for(unsigned int i = 0; i<j.size(); i++){
	j[i] = A.j[i];
	a[i] = A.a[i];
      }
    }
  
  
    Vecteur addmatmul(Vecteur const& x, Vecteur & Ax) const
    {
        for(int i=0; i< n; ++i)
            for( int k=ip[i]; k< ip[i+1]; ++k)
	      Ax(i) += a[k]*x(j[k]);
        return Ax;
    }
    double * pij(int I,int J)
    {   // dichotimie
        int k0=ip[I],k1=ip[I+1]-1, k ;
        // recheche J dans k = J k = ip[I], ip[I+1]=1 :
        while (k0<=k1)
        {
            k=(k0+k1)/2;
            if( j[k]==J) return &(a[k]);
            else if (j[k] > J) k1=k-1;
            else  k0=k+1;
        }
        return 0;
    }
    double & operator()(int i, int j)
    { double *paij=pij(i,j); assert(paij); return *paij;}
};





std::ostream &   operator<<(std::ostream & f, const std::vector<double> & b);
std::ostream &   operator<<(std::ostream & f, const MatCSR & A);



