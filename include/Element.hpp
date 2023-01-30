#ifndef Element_HPP
#define Element_HPP

#include"Vertex.hpp"

class Element{

public:
  Element();
  Element(Element const& v);
  ~Element();

   /*visualization methods*/
  void print_element() const;

  /*getters*/
  double get_vertices(int i) const;
  double get_label() const;
 
 


private:
  int vertices[3];
  int label;



}

#endif
