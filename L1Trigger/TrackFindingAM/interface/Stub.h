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

  void setPhi(const double & phi) { phi_ = phi; }
  void setR(const double & R) { R_ = R; }
  void setZ(const double & z) { z_ = z; }
  double phi() const { return phi_; }
  double R() const { return R_; }
  double z() const { return z_; }
  int layer() const { return layer_; }
  float strip() const { return strip_; }

 private:
  double phi_;
  double R_;
  double z_;
  int layer_;
  float strip_;
};


#endif //REMOTEPROJECTS_STUB_H
