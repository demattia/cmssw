#ifndef _CMSPATTERNLAYER_H_
#define _CMSPATTERNLAYER_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <bitset>
#include "PatternLayer.h"
#include <stdexcept>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/map.hpp>

#ifdef IPNL_USE_CUDA
#include "gpu_struct.h"
#endif

using namespace std;


/**
   \brief First version of a CMS pattern structure (16 bits are used)
   Also contains all the detector geometry :
     - ids of layers
     - number of ladders per layer
     - number of modules per ladder
     - number of strips per segment
     - convertion from root file module ID to program module ID
     - convertion from root file ladder ID to program ladder ID
   Please keep in mind that the value 15 for the ladder is used for the fake stubs.
   The superstrip position value is stored in Gray code, this allows to manage more interesting situations with the DC bits.
**/

class CMSPatternLayer : public PatternLayer{
 private:

  short binaryToGray(short num);
  short grayToBinary(short gray);

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
    {
      ar << boost::serialization::base_object<PatternLayer>(*this);
    }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
    {
      ar >> boost::serialization::base_object<PatternLayer>(*this);
      if(version<2){
	cout<<"This release does not support old PBK files. Please use up-to-date files or switch to an older release."<<endl;
	exit(-1);
      }
      if(version<1){
	//Convert the strip value to gray code
	char current_val = getStripCode();
	bool even = (current_val%2)==0;
	bits &= ~(STRIP_MASK<<STRIP_START_BIT);
	bits |= (binaryToGray(current_val)&STRIP_MASK)<<STRIP_START_BIT;
	if(!even){
	  // Update the DC bits (change the first DC bits only if the strip is even)
	  if(dc_bits[0]<2)
	    dc_bits[0]=1-dc_bits[0];
	}
      }
    }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

 public:

  static short MOD_START_BIT;
  static short PHI_START_BIT;
  static short STRIP_START_BIT;
  static short SEG_START_BIT;

  static short MOD_MASK;
  static short PHI_MASK;
  static short STRIP_MASK;
  static short SEG_MASK;

  static short OUTER_LAYER_SEG_DIVIDE;//Simplification factor on outer barrel layers segments (1:we use segment, 2:all values to 0)
  static short INNER_LAYER_SEG_DIVIDE;//Simplification factor on inner barrel layers segments (1:we use segment, 2:all values to 0)

  CMSPatternLayer();
  CMSPatternLayer* clone();
  vector<SuperStrip*> getSuperStrip(int l, Detector& d);
  void getSuperStripCuda(int l, const vector<int>& ladd, const map<int, vector<int> >& modules, int layerID, unsigned int* v);
  
  /**
     \brief Set the values in the patternLayer
     \param m The module Z position (0 to 13 for modules 23 to 47)
     \param phi The phi position of the module in the sector (0 to 7)
     \param strip The super strip number
     \param seg The segment in the module (0 or 1)
  **/
  void setValues(short m, short phi, short strip, short seg);

  /**
     \brief Compute the superstrip value from the stub informations
     \param layerID The ID of the stub's layer
     \param module The value of the stub's module in the trigger tower
     \param phi The value of the stub's ladder in the trigger tower
     \param strip The value of the stub's strip
     \param seg The value of the stub's segment
     \param sstripSize The size of a superstrip in numbner of strips
     \param isPS True if we are on a P/S module
     \param fake True if you are building a fake superstrip
   **/
  void computeSuperstrip(short layerID, short module, short phi, short strip, short seg, int sstripSize, bool isPS, bool fake=0);

  /**
     \brief Returns a string representation of the PatternLayer
     \return A string describing the PatternLayer
  **/
  string toString();

  /**
     \brief Returns a string representation of the PatternLayer, using binary values
     \return A string describing the PatternLayer
  **/
  string toStringBinary();
  /**
     \brief Returns a string representation of the PatternLayer, using the encoding needed for a AM05 chip
     \param tagLayer If true, set the bit 6 to 1 (used to distinguish layers in case of 2 layers on 1 bus)
     \return A string describing the PatternLayer
  **/
  string toAM05Format(bool tagLayer=false);
  /**
     \brief Returns the module's Z position
     \return The module's Z position
  **/
  short getModule();
  /**
     \brief Returns the ladder phi position
     \return The ladder's phi position
  **/
  short getPhi();
  /**
     \brief Returns the Super strip position
     \return The position of the super strip in the segment
  **/
  short getStrip();
  /**
     \brief Returns the Super strip encoded value
     \return The gray encoded value of the position of the super strip in the segment
  **/
  short getStripCode();
  /**
     \brief Returns the position of the segment in the module
     \return The segment's position in the module (0 or 1)
  **/
  short getSegment();

  /**
     \brief Returns the list of layers IDs in the detector
     \return a vector containing the layer IDs
  **/
  static vector<int> getLayerIDs();

  /**
     \brief Returns the number of strips in a segment
     \return an int value
  **/
  static int getNbStripsInSegment();

  /**
     \brief Get the ID of the ladder from the ladder ID in the muon file. Used to change the IDs between the root file and the simulation program (if needed).
     \brief After a call to this method, ladders numbering must start at 0
     \param layerID The layer ID (TIB : 5,6,7 - TOB : 8,9,10 - TEC : 11,12,13,14,15,16,17 and 18,19,20,21,22,23,24)
     \param ladderID The ladder ID in the muon file
     \return The ID to use in the program
  **/
  static int getLadderCode(int layerID, int ladderID);

 /**
     \brief Get the ID of the module from the module ID in the muon file. Used to change the IDs between the root file and the simulation program (if needed).
     \brief This method allows to change the module resolution (you can divide the ID by 2) or the numbering.
     \brief After a call to this method, modules numbering must start at 0
     \param layerID The layer ID (TIB : 5,6,7 - TOB : 8,9,10 - TEC : 11,12,13,14,15,16,17 and 18,19,20,21,22,23,24)
     \param moduleID The module ID in the muon file
     \return The ID to use in the program
  **/
  static int getModuleCode(int layerID, int moduleID);

  /**
     \brief Get the code of the segment in the patternLayer from the segment ID in the muon file
     \brief Segment ID must be 0 or 1
     \param layerID The layer ID (TIB : 5,6,7 - TOB : 8,9,10 - TEC : 11,12,13,14,15,16,17 and 18,19,20,21,22,23,24)
     \param ladderID The ladder ID in the muon file
     \param segmentID The segment ID in the muon file
     \return The number to use in the simulation program
  **/
  static int getSegmentCode(int layerID, int ladderID, int segmentID);

  /**
     \brief Gives the number of ladders in the detector for a given layer
     \param layerID The layer ID
     \return The number of ladders on the layer layerID
  **/
  static int getNbLadders(int layerID);

  /**
     \brief Gives the number of modules in the detector for a given layer/ladder
     \param layerID The layer ID
     \param ladderID The ladder ID
     \return The number of modules on the layer layerID / ladder ladderID
  **/
  static int getNbModules(int layerID, int ladderID);  

  /**
     \brief Convertion from CMSSW layer numbering (one id per layer/disk) to PRBF2 numbering (id depends on the layer AND the module type)
     \param cms_layer The CMSSW layer ID (5 to 10 for barrel, 11 to 15 for first TEC and 18 to 22 for second TEC)
     \param isPS True if the module is PS, False if 2S
     \return The PRBF2 layer ID
  **/
  static int cmssw_layer_to_prbf2_layer(int cms_layer, bool isPS);

 /**
     \brief Check if the PatternLayer is a fake one (used on layers not crossed by the track)
     \return True if the PatternLayer is a placeholder
  **/  
  bool isFake();

  vector<int> getHDSuperstrips();  

/**
     \brief Returns a map containing the valid ETA range for each layer
     \return For each layerID, gives the minimum and maximum ETA values
  **/
  static map<int, pair<float,float> > getLayerDefInEta();
};
BOOST_CLASS_VERSION(CMSPatternLayer, 2)
#endif
