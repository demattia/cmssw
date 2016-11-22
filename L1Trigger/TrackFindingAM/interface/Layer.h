#ifndef _LAYER_H_
#define _LAYER_H_

#include "Ladder.h"
#include <vector>

/**
   \brief Represents a Layer (vector of Ladders)
**/
class Layer{

 private:
  vector<Ladder*> ladders;
  bool isBarrelLayer;
  
 public:
  /**
     \brief Constructor
     \param nbLad Number of ladders on this layer
     \param nbMod Number of modules on each ladder
     \param nbSeg Number of segments on each module
     \param segmentSize Number of strips in a segment
     \param sstripSize Number of strips in a super strip
     \param isBarrel True if it's a barrel layer
  **/
  Layer(int nbLad, int nbMod, int nbSeg, int segmentSize, int sstripSize, bool isBarrel);
  /**
     \brief Destructor
  **/
  ~Layer();
  /**
     \brief Get one of the ladders
     \param pos The ladder ID starting from 0.
     \return A pointer on the Ladder object (not a copy), NULL if not found.
  **/
  Ladder* getLadder(int pos);
  /**
     \brief Desactivates all the superstrips in the layer.
  **/
  void clear();

  bool isBarrel();

};
#endif
