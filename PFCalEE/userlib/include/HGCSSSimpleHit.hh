#ifndef _hgcsssimplehit_hh_
#define _hgcsssimplehit_hh_

#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <iostream>


class HGCSSSimpleHit{
  
public:
  
  HGCSSSimpleHit():
    e_(0),
    l_(0),
    x_(0),
    y_(0),
    z_(0)
  {};
    

  ~HGCSSSimpleHit(){};
    
  inline double e(){
    return e_;
  };
  inline double x(){
    return x_;
  };
  inline double y(){
    return y_;
  };
  inline double z(){
    return z_;
  };
  inline unsigned l(){
    return l_;
  };

  inline void setE(const double & val){
    e_ = val;
  };
  inline void setx(const double & val){
    x_ = val;
  };
  inline void sety(const double & val){
    y_ = val;
  };
  inline void setz(const double & val){
    z_ = val;
  };
  inline void setLayer(const unsigned val){
    l_ = val;
  };

 private:
  double e_;
  unsigned l_;
  double x_;
  double y_;
  double z_;


  ClassDef(HGCSSSimpleHit,1);


};//class

typedef std::vector<HGCSSSimpleHit> HGCSSSimpleHitVec;

#endif
