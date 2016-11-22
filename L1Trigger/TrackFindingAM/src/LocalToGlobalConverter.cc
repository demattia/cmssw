#include "../interface/LocalToGlobalConverter.h"

LocalToGlobalConverter::LocalToGlobalConverter(const Sector& sectorDefinition, string geometryFile){
  string line;
  ifstream myfile (geometryFile.c_str());
  int layer = -1;
  int ladder = -1;
  int local_ladder = -1;
  int module = -1;
  int local_module = -1;
  float coef_value = -1;
  stringstream val;

  module_pos = NULL;

  //Which side of the tracker are we on?
  if(sectorDefinition.getOfficialID()<24)
    tracker_side=true;
  else
    tracker_side=false;

  //Parse the geometry file containing the cartesian coordinates of all modules
  if (myfile.is_open()){
        
    //Memory allocation of the structure
    module_pos = new vector<float>**[16];
    for(int i=0;i<16;i++){//we can have up to 16 layers
      module_pos[i] = new vector<float>*[16];
      for(int j=0;j<16;j++){//16 ladders
	module_pos[i][j] = new vector<float>[32];
	for(int k=0;k<32;k++){//32 modules
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	}
      }
    }
    
    while ( myfile.good() ){
      getline (myfile,line);
      if(line.length()>0 && line.find("#")!=0){
	stringstream ss(line);
	std::string item;
	vector<string> items;
	while (getline(ss, item, '/')) {
	  std::string::iterator end_pos = std::remove(item.begin(), item.end(), ' ');
	  item.erase(end_pos, item.end());
	  items.push_back(item);
	}
	if(items.size()==11){
	  val.clear();
	  val.str(items[0]);
	  val >> layer;
	  val.clear();
	  val.str(items[1]);
	  val >> ladder;
	  ladder = ladder-1;//numbering in file starts with 1 -> corrects to 0
	  val.clear();
	  val.str(items[2]);
	  val >> module;
	  local_ladder = sectorDefinition.getLadderCode(layer,ladder);
	  if(local_ladder==-1)//not in the sector
	    continue;
	  local_module = sectorDefinition.getModuleCode(layer,ladder,module);
	  if(local_module==-1)//not in the sector
	    continue;
	  bool isPSModule = false;
	  if((layer>=5 && layer<=7) || (layer>10 && ladder<=8)){
	    isPSModule=true;
	  }
	  int prbf2_layer = CMSPatternLayer::cmssw_layer_to_prbf2_layer(layer,isPSModule);

	  val.clear();
	  val.str(items[8]);
	  val >> coef_value;
	  module_pos[prbf2_layer][local_ladder][local_module][0]=coef_value;
	  val.clear();
	  val.str(items[9]);
	  val >> coef_value;
	  module_pos[prbf2_layer][local_ladder][local_module][1]=coef_value;
	  val.clear();
	  val.str(items[10]);
	  val >> coef_value;
	  module_pos[prbf2_layer][local_ladder][local_module][2]=coef_value;

	  //Computes the angle between X axis and the module center
	  float PhiMod = atan2(module_pos[prbf2_layer][local_ladder][local_module][1],module_pos[prbf2_layer][local_ladder][local_module][0]);//atan2(y,x)
	  float fModuleWidth,fModuleHeight,nModuleStrips,nModuleSegments;
	  if (isPSModule){
	    //PS
	    fModuleWidth 	= 4.48144;
	    fModuleHeight 	= 9.59;
	    nModuleStrips	= 959.0;
	    nModuleSegments	= 31.0;
	  }
	  else{
	    //2S
	    fModuleWidth 	= 5.025;
	    fModuleHeight 	= 9.135;
	    nModuleStrips	= 1015.0;
	    nModuleSegments	= 1.0;
	  }
	  float fStripPitch = fModuleHeight / nModuleStrips;
	  float fSegmentPitch = fModuleWidth / nModuleSegments;

	  if(layer<11){
	    module_pos[prbf2_layer][local_ladder][local_module][3]=fStripPitch * sin(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][4]=-fStripPitch * cos(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][5]=0;
	    module_pos[prbf2_layer][local_ladder][local_module][6]=0;
	    module_pos[prbf2_layer][local_ladder][local_module][7]=0;
	    module_pos[prbf2_layer][local_ladder][local_module][8]=-fSegmentPitch;
	  }
	  else{
	    module_pos[prbf2_layer][local_ladder][local_module][3]=fStripPitch * sin(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][4]=-fStripPitch * cos(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][5]=0;
	    module_pos[prbf2_layer][local_ladder][local_module][6]=-fSegmentPitch * cos(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][7]=-fSegmentPitch * sin(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][8]=0;
	  }
	}
      }
    }
    myfile.close();
  }
  else{
    cout << "Can not find file "<<geometryFile<<" to load the modules position lookup table!"<<endl;
  }
}

LocalToGlobalConverter::~LocalToGlobalConverter(){ 
  if(module_pos!=NULL){
    for(int i=0;i<16;i++){//we can have up to 16 layers
      for(int j=0;j<16;j++){//16 ladders
	for(int k=0;k<32;k++){//32 modules
	  module_pos[i][j][k].clear();
	}
	delete [] module_pos[i][j];
      }
      delete [] module_pos[i];
    }
    delete [] module_pos;
  }
}

vector<float> LocalToGlobalConverter::toGlobal(int layer, int ladder, int module, int segment, float strip) const throw (std::runtime_error){

  if(module_pos==NULL){
    throw std::runtime_error("Modules position lookup table not found");
  }

  bool isPS = (layer<8);
  float relatStrip=0;
  float relatSeg=0;
  float X;
  float Y;
  
  bool isBarrel = ((layer&0x7)<3);

  if(isPS){
    relatStrip = strip-(959.0/2.0);
    relatSeg   = segment-(31.0/2.0);
  }
  else{
    relatStrip = strip-(1015.0/2.0);
    relatSeg   = segment-(1.0/2.0);
  }

  vector<float> res;
  vector<float> positions = module_pos[layer][ladder][module];

  if(!tracker_side && !isBarrel){
    X = positions[0] - relatStrip*positions[3] + relatSeg*positions[6];
    Y = positions[1] - relatStrip*positions[4] + relatSeg*positions[7];
  }
  else{
    X = positions[0] + relatStrip*positions[3] + relatSeg*positions[6];
    Y = positions[1] + relatStrip*positions[4] + relatSeg*positions[7];
  }
  res.push_back(X);//X
  res.push_back(Y);//Y
  res.push_back(positions[2] + relatStrip*positions[5] + relatSeg*positions[8]);//Z

  return res;
}
