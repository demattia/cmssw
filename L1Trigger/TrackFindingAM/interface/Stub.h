//
// Created by Marco De Mattia on 7/30/15.
//

#ifndef REMOTEPROJECTS_STUB_H
#define REMOTEPROJECTS_STUB_H


class Stub
{
 public:
  Stub(const double & phi, const double & R, const double & z, const int layer, const int strip) :
      phi_(phi), R_(R), z_(z), layer_(layer), strip_(strip)
  {}

  Stub(const double & phi, const double & R, const double & z, const int layer, const int strip,
       const unsigned int stubRef) :
      phi_(phi), R_(R), z_(z), layer_(layer), strip_(strip), stubRef_(stubRef)
  {}

  void setPhi(const double & phi) { phi_ = phi; }
  void setR(const double & R) { R_ = R; }
  void setZ(const double & z) { z_ = z; }
  void setStubRef(const unsigned int stubRef) { stubRef_ = stubRef; }
  double phi() const { return phi_; }
  double R() const { return R_; }
  double z() const { return z_; }
  int layer() const { return layer_; }
  float strip() const { return strip_; }
  unsigned int stubRef() const { return stubRef_; }

 private:
  double phi_;
  double R_;
  double z_;
  int layer_;
  float strip_;
  // Only needed for the CombinationBuilder
  unsigned int stubRef_;
};


#endif //REMOTEPROJECTS_STUB_H
