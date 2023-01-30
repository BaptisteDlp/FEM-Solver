#ifndef MESHS_HPP
#define MESHS_HPP


#include<vector>
#include<cassert>
#include<utility>
#include"Vecteur.hpp"





class Mesh{

public:
  Mesh(std::string file_name, bool offset = false);
  ~Mesh();

  void uniform_mesh(int nbv_h, int nbv_v, double x0, double xN, double y0, double yN, std::string fileName = "Th.msh");

  void print_vertices() const;
  void print_elements() const;
  void print_edges() const;
  void print_vertices_boundary() const;

  /*visualization methods*/
  void export_mesh_data() const;
  void mesh_visu(bool mesh = true, bool solution = false) const;
  void export_solution_data(Vecteur const& solution) const;
 

  /*getters*/
  int get_nbv() const;
  int get_nbt() const;
  int get_nbe() const;
  int get_nbeb() const;
  int get_nbvBound() const;
  double get_x(int index) const;
  double get_y(int index) const;
  int get_vertex(int index, int v) const;
  int get_edge(int index, int v) const;
  int get_vertex_bound(int index) const;
  /*setters*/
  void set_resolved();
  void unset_resolved();

  friend void  build_edges(double **tabv, int nbv, int nbt, int *nbe, int *nbeb, int *nbvBound, int **tabt, std::vector<std::pair<int,int>> &tabe, std::vector<std::pair<int,int>> &tabeb, std::vector<int> &tabvBound);
  
  

private:
  std::string mesh_name;
  std::string mesh_data;
  std::string mesh_gnu;

  std::string solution_data;
  std::string solution_gnu;

  bool resolved;
  
  double **tabv;
  int **tabt;
  std::vector<std::pair<int,int>> tabe;
  std::vector<std::pair<int,int>> tabeb;
  std::vector<int> tabvBound;
  int nbv;
  int nbt;
  int nbe;
  int nbeb;
  int nbvBound;
  

};



#endif
