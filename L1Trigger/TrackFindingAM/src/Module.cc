#include "../interface/Module.h"

Module::Module(int nbSeg, int segmentSize, int sstripSize){
  for(int i=0;i<nbSeg;i++){
    Segment* s = new Segment(segmentSize, sstripSize);
    segments.push_back(s);
  }
}

Module::~Module(){
  for(unsigned int i=0;i<segments.size();i++){
    delete segments[i];
  }
  segments.clear();
}

Segment* Module::getSegment(int n){
  if(n>-1 && n<(int)segments.size())
    return segments[n];
  return NULL;
}

void Module::clear(){
  for(unsigned int i=0;i<segments.size();i++){
    segments[i]->clear();
  }
}
