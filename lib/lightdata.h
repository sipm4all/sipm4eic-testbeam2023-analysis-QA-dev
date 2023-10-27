#pragma once

namespace sipm4eic {

class lightdata {

 public:

  unsigned char device = 0;
  unsigned char index = 0;
  unsigned char coarse = 0;
  unsigned char fine = 0;

  int chip() { return index / 32; };
  int eoch() { return index % 64; };
  
 lightdata(unsigned char _device, unsigned char _index, unsigned char _coarse, unsigned char _fine) :
  device(_device),
    index(_index),
    coarse(_coarse),
    fine(_fine) {
 }
    
    
  
};

} /** namespace sipm4eic **/
