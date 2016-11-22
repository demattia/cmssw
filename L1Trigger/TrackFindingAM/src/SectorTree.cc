#include "../interface/SectorTree.h"

map<string, int> SectorTree::superstripSize_lut;
string SectorTree::ss_size_filename="./ss_size.txt";

SectorTree::SectorTree(){
  srand ( time(NULL) );
  mapNeedsUpdate=true;
}

SectorTree::SectorTree(const SectorTree& st){
  srand ( time(NULL) );
  mapNeedsUpdate=true;
  //copy the superstrip sizes
  for(map<string, int>::iterator it=st.superstripSize_lut.begin();it!=superstripSize_lut.end();it++){
    this->superstripSize_lut[it->first]=it->second;
  }
}

SectorTree::~SectorTree(){
  for(unsigned int i=0;i<sector_list.size();i++){
    delete sector_list[i];
  }
}

Sector* SectorTree::getSector(vector<int> ladders, vector<int> modules){

  // check that the multimap is populated
  if(mapNeedsUpdate)
    updateSectorMap();

  pair<multimap<string,Sector*>::iterator,multimap<string,Sector*>::iterator> ret;
  multimap<string,Sector*>::iterator first;

  ostringstream oss;
  oss<<std::setfill('0');
  for(unsigned int j=0;j<ladders.size();j++){
    oss<<setw(2)<<ladders[j];
  }

  ret = sectors.equal_range(oss.str());

  if(ret.first==ret.second){//NOT FOUND
    return NULL;
  }

  first=ret.first;  
  while(first!=ret.second){
    Sector* test = first->second;
    
    bool found=false;
    for(unsigned int i=0;i<ladders.size();i++){
      vector<int> mods = test->getModules(i,ladders[i]);
      found = false;
      for(unsigned int j=0;j<mods.size();j++){
	if(modules[i]==mods[j]){
	  found=true;
	  break;
	}
      }
      if(!found)
	break;
    }
    if(found)//this one is ok-> we take it
      return test;
    first++;
  }
  //none of the selected sectors where ok for modules...
  return NULL;
}

Sector* SectorTree::getSector(const Hit& h){
  for(unsigned int i=0;i<sector_list.size();i++){
    if(sector_list[i]->contains(h))
      return sector_list[i];
  }
  return NULL;
}

void SectorTree::addSector(Sector s){
  Sector* ns = new Sector(s);
  sector_list.push_back(ns);
  mapNeedsUpdate = true;
}

void SectorTree::updateSectorMap(){
  sectors.clear();
  for(unsigned int j=0;j<sector_list.size();j++){
    Sector* ns = sector_list[j];
    vector<string> keys = ns->getKeys();
    for(unsigned int i=0;i<keys.size();i++){
      sectors.insert(pair<string, Sector*>(keys[i],ns));
    }
  }
  mapNeedsUpdate=false;
}

vector<Sector*> SectorTree::getAllSectors(){
  return sector_list;
}

int SectorTree::getLDPatternNumber(){
  int nb = 0;
  for(unsigned int i=0;i<sector_list.size();i++){
    nb+=sector_list[i]->getLDPatternNumber();
  }
  return nb;
}

int SectorTree::getFDPatternNumber(){
  int nb = 0;
  for(unsigned int i=0;i<sector_list.size();i++){
    nb+=sector_list[i]->getFDPatternNumber();
  }
  return nb;  
}

void SectorTree::computeAdaptativePatterns(map<int, int> r){
  for(unsigned int i=0;i<sector_list.size();i++){
    Sector* s=sector_list[i];
    s->computeAdaptativePatterns(r);
  }
}

void SectorTree::link(Detector& d){
  for(unsigned int i=0;i<sector_list.size();i++){
    sector_list[i]->link(d);
  }
}

#ifdef IPNL_USE_CUDA
void SectorTree::linkCuda(patternBank* p, deviceDetector* d){
  for(unsigned int i=0;i<sector_list.size();i++){
    sector_list[i]->linkCuda(p,d);
  }
}
#endif

bool comparePatternOrder(GradedPattern* p1, GradedPattern* p2){
  return p1->getOrderInChip()<p2->getOrderInChip();
}

vector<Sector*> SectorTree::getActivePatternsPerSector(int active_threshold, unsigned int max_nb_roads){
  vector<Sector*> list;
  for(unsigned int i=0;i<sector_list.size();i++){
    Sector* copy = new Sector(*sector_list[i]);
    vector<GradedPattern*> active_patterns = sector_list[i]->getActivePatterns(active_threshold);
    sort(active_patterns.begin(),active_patterns.end(),comparePatternOrder);//order the roads by their chip's address
    for(unsigned int j=0;j<active_patterns.size();j++){
      if(j<max_nb_roads)
	copy->getPatternTree()->addPattern(active_patterns[j], NULL, active_patterns[j]->getAveragePt());
      delete active_patterns[j];
    }
    list.push_back(copy);
  }
  return list;
}

vector<Sector*> SectorTree::getActivePatternsPerSectorUsingMissingHit(int max_nb_missing_hit, int active_threshold, unsigned int max_nb_roads){
  vector<Sector*> list;
  for(unsigned int i=0;i<sector_list.size();i++){
    Sector* copy = new Sector(*sector_list[i]);
    vector<GradedPattern*> active_patterns = sector_list[i]->getActivePatternsUsingMissingHit(max_nb_missing_hit, active_threshold);
    sort(active_patterns.begin(),active_patterns.end(),comparePatternOrder);//order the roads by their chip's address
    for(unsigned int j=0;j<active_patterns.size();j++){
      if(j<max_nb_roads)
	copy->getPatternTree()->addPattern(active_patterns[j], NULL, active_patterns[j]->getAveragePt());
      delete active_patterns[j];
    }
    list.push_back(copy);
  }
  return list;
}

map< string, int > SectorTree::loadSStripSizeLUT(string name){
  string line;
  ifstream myfile (name.c_str());
  map< string, int > size_lut;
  if (myfile.is_open()){
    while ( myfile.good() ){
      getline (myfile,line);
      if(line.length()>0 && line.find("#")!=0){
	stringstream ss(line);
	std::string item;
	vector<string> items;
	while (getline(ss, item, ' ')) {
	  std::string::iterator end_pos = std::remove(item.begin(), item.end(), ' ');
	  item.erase(end_pos, item.end());
	  items.push_back(item);
	}
	if(items.size()==2){
	  istringstream buffer2(items[1]);
	  int size;
	  buffer2 >> size;
	  size_lut[items[0]]=size;
	}
      }
    }
    myfile.close();
  }
  else{
    cout << "Can not find file "<<name<<" to load the superstrip size lookup table!"<<endl;
    exit(-1);
  }
  return size_lut;
}

int SectorTree::getSuperstripSize(int layer_id, int ladder_id){
  if(superstripSize_lut.size()==0)
    superstripSize_lut = loadSStripSizeLUT(ss_size_filename);
  if(layer_id<11){ // barrel
    ostringstream oss;
    oss<<std::setfill('0');
    oss<<setw(2)<<(int)layer_id;
    return superstripSize_lut[oss.str()];
  }
  else{//endcap
    ostringstream oss;
    oss<<std::setfill('0');
    oss<<setw(2)<<(int)layer_id;
    oss<<setw(2)<<(int)ladder_id;
    return superstripSize_lut[oss.str()];
  }
}

map<string, int> SectorTree::getSuperstripSize_lut(){
  if(superstripSize_lut.size()==0)
    superstripSize_lut = loadSStripSizeLUT(ss_size_filename);
  return superstripSize_lut;
}

bool SectorTree::hasSameSuperstripSizes(const SectorTree& st){
  //all of this is in st
  for(map<string, int>::iterator it=superstripSize_lut.begin();it!=superstripSize_lut.end();it++){
    if(it->second!=st.superstripSize_lut[it->first])
      return false;
  }
  //all of st is in this
  for(map<string, int>::iterator it=st.superstripSize_lut.begin();it!=st.superstripSize_lut.end();it++){
    if(it->second!=superstripSize_lut[it->first])
      return false;
  }
  return true;
}

void SectorTree::displaySuperstripSizes(){
  for(map<string, int>::iterator it=superstripSize_lut.begin();it!=superstripSize_lut.end();it++){
    cout<<"\t"<<it->first<<" : "<<it->second<<endl;
  }
}

int SectorTree::getNbSectors(){
  return sector_list.size();
}

void SectorTree::setSuperstripSizeFile(string fileName){
  ss_size_filename=fileName;
}
