#ifndef Vertex_HPP
#define Vertex_HPP

#include"R2.hpp"

class Vertex{

public:
  Vertex();
  Vertex(Vertex const& v);
  ~Vertex();

   /*visualization methods*/
  void print_vertex() const;

  /*getters*/
  double get_coord() const;
  double get_label() const;
 
 


private:
  R2 coord;
  int label;



}

#endif
