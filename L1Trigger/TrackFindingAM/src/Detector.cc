#include "../interface/Detector.h"

#include "../interface/SectorTree.h"

Detector::Detector(){
  dump=NULL;
  verbose = false;
}

void Detector::addLayer(int lNum, int nbLad, int nbMod, int nbSeg, int segmentSize, int sstripSize, bool barrel){
  if(dump==NULL)
      dump = new SuperStrip(sstripSize);
  Layer* l = new Layer(nbLad, nbMod, nbSeg, segmentSize, sstripSize, barrel);
  layerNumber.push_back(lNum);
  superStripSizes.push_back(sstripSize);
  layers.push_back(l);
}

void Detector::setSectorMaps(map<string,int> lm, map<string,int> mm){
  ladderMap = lm;
  moduleMap = mm;
}

Detector::~Detector(){
  if(dump!=NULL)
    delete dump;
  for(unsigned int i=0;i<layers.size();i++){
    delete layers[i];
  }
  layers.clear();
  layerNumber.clear();
  ladderMap.clear();
  moduleMap.clear();
}

int Detector::getLayerPosition(int pos){
  for(unsigned int i=0;i<layerNumber.size();i++){
    if(layerNumber[i]==pos)
      return i;
  }
  return -1;
}

Layer* Detector::getLayerFromAbsolutePosition(int pos){
  if(pos>-1 && pos<(int)layers.size())
    return layers[pos];
  return NULL;
}

Layer* Detector::getLayer(int pos){
  int localPosition = getLayerPosition(pos);
  if(localPosition!=-1)
    return layers[localPosition];
  return NULL;
}

void Detector::clear(){
  for(unsigned int i=0;i<layers.size();i++){
    layers[i]->clear();
  }
}

void Detector::receiveHit(const Hit& h){
  if(verbose)
    cout<<"#STUB : "<<h<<endl;
  int l = getLayerPosition(h.getLayer());
  if(l!=-1){
    Layer* la = getLayerFromAbsolutePosition(l);
    if(la!=NULL){
      try{
	CMSPatternLayer pat;

	bool isPSModule = false;
	if((h.getLayer()>=5 && h.getLayer()<=7) || (h.getLayer()>10 && h.getLadder()<=8)){
	  isPSModule=true;
	}
	
	ostringstream oss;
	oss<<std::setfill('0');
	oss<<setw(2)<<(int)h.getLayer();
	oss<<setw(2)<<(int)h.getLadder();
	string lad = oss.str();
	oss<<setw(2)<<(int)h.getModule();

	pat.computeSuperstrip(h.getLayer(), moduleMap[oss.str()], ladderMap[lad], h.getStripNumber(), h.getSegment(), SectorTree::getSuperstripSize(h.getLayer(),h.getLadder()), isPSModule);
	if(verbose){
	  cout<<"#SUPERSTRIP : "<<(int)h.getLayer()<<" "<<hex<<"0x"<<std::setfill ('0') << std::setw (4)<<pat.getIntValue()<<dec<<endl;
	  cout<<endl;
	}
	SuperStrip* s = la->getLadder(pat.getPhi())->getModule(la->isBarrel()?0:pat.getModule())->getSegment(la->isBarrel()?pat.getModule()*2+pat.getSegment():pat.getSegment())->getSuperStripFromIndex(pat.getStrip());

	if(s==NULL)
	  cout<<"ERROR : Cannot find superStrip corresponding to the following hit : "<<h<<endl;
	else
	  s->touch(&h);
      }
      catch (const std::out_of_range& oor) {
	std::cerr << "The following point cannot be mapped to a superstrip : "<<endl;
	std::cerr << h << endl;
      }
    }
  }
  else{
    //int tmp_layer = (int)h.getLayer();
    //cout<<"no layer "<<tmp_layer<<endl;
  }
}

int Detector::getNbLayers(){
  return (int)layers.size();
}

SuperStrip* Detector::getDump(){
  return dump;
}

void Detector::setVerboseMode(bool m){
  verbose = m;
}
