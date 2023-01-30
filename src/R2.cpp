#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<cstring>
#include<cmath>
#include"R2.hpp"


using namespace std;

double R2::dist(R2 const& v) const {

  return sqrt(pow(x-v.get_x(),2)+pow(y-v.get_y(),2));

}

ostream& operator<<(std::ostream &os, R2 const& v){

  os<<"("<<v.get_x()<<","<<v.get_y()<<")"<<endl;
}
