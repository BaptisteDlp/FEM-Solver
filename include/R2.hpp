#ifndef R2_HPP
#define R2_HPP




class R2{

public:
  R2() : x(0.), y(0.){};
  R2(double x1, double x2) : x(x1), y(x2){};
  R2(R2 const& v) : x(v.get_x()), y(v.get_y()){};
  ~R2();


  double dist(R2 const& v) const;
  
  /*visualization methods*/
  void print_point() const;

  /*getters*/
  double get_x() const;
  double get_y() const;
 
  /*setters*/
  double set_x();
  double set_y();
  
  

private:
  double x;
  double y;

  

};


std::ostream& operator<<(std::ostream &os, R2 const& v);
		      
#endif
