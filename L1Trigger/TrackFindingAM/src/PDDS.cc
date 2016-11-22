//
// Created by Denis Rathjens on 11/17/16.
//

#include "L1Trigger/TrackFindingAM/interface/PDDS.h"

PDDS::PDDS(){}

unsigned int PDDS::IndexInterpreter(unsigned int Index_){ //interprets the index of Marco's 6/6 combinations for easy use

/*Special combinations printer
  std::cout<<"[5,0],[6,0],[7,0],[8,0]: "<<combinationIndex({5,6,7,8},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[12,0],[13,0],[14,0],[15,0]: "<<combinationIndex({12,13,14,15},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[11,0],[12,0],[13,0],[14,0]: "<<combinationIndex({11,12,13,14},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[9,0]: "<<combinationIndex({5,6,8,9},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[11,0],[13,0],[14,0],[15,0]: "<<combinationIndex({11,13,14,15},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[11,0],[12,0],[13,0],[15,0]: "<<combinationIndex({11,12,13,15},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[8,0],[9,0]: "<<combinationIndex({6,7,8,9},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[9,0]: "<<combinationIndex({5,6,7,9},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[6,0],[11,70],[12,70],[13,110],[14,110]: "<<combinationIndex({6,11,12,13,14},{0.,70.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[12,70],[13,110],[14,110]: "<<combinationIndex({5,6,11,12,13,14},{0.,0.,70.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,0],[12,0],[13,0],[14,70]: "<<combinationIndex({5,6,11,12,13,14},{0.,0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[6,0],[11,0],[12,0],[13,0],[14,70]: "<<combinationIndex({6,11,12,13,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[11,0],[12,0],[14,0],[15,0]: "<<combinationIndex({11,12,14,15},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[11,110]: "<<combinationIndex({5,6,7,11},{0.,0.,0.,110.},14)<<std::endl;
  std::cout<<"[5,0],[7,0],[8,0],[9,0]: "<<combinationIndex({5,7,8,9},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[11,110]: "<<combinationIndex({5,6,8,11},{0.,0.,0.,110.},14)<<std::endl;
  std::cout<<"[5,0],[11,70],[12,70],[13,110],[14,110]: "<<combinationIndex({5,11,12,13,14},{0.,70.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[8,0],[11,110]: "<<combinationIndex({6,7,8,11},{0.,0.,0.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[13,110],[14,110]: "<<combinationIndex({5,11,12,13,14},{0.,0.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[12,110]: "<<combinationIndex({5,6,7,12},{0.,0.,0.,110.},14)<<std::endl;
  std::cout<<"[5,0],[7,0],[8,0],[11,110]: "<<combinationIndex({5,7,8,11},{0.,0.,0.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,0],[12,0],[13,0]: "<<combinationIndex({5,6,11,12,13},{0.,0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[12,70],[14,110]: "<<combinationIndex({5,6,11,12,14},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[11,70],[12,110],[13,110]: "<<combinationIndex({5,6,8,11,12,13},{0.,0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[12,0],[13,0],[14,70]: "<<combinationIndex({5,6,12,13,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[11,0],[12,0],[13,70],[14,70],[15,110]: "<<combinationIndex({11,12,13,14,15},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[13,110]: "<<combinationIndex({5,6,11,13},{0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,0],[13,0],[14,70]: "<<combinationIndex({5,6,11,13,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[11,70]: "<<combinationIndex({5,6,7,11},{0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[6,0],[11,70],[12,70],[13,110]: "<<combinationIndex({6,11,12,13},{0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[12,110]: "<<combinationIndex({5,6,11,12},{0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[12,70],[13,110]: "<<combinationIndex({5,6,12,13},{0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,110],[12,110]: "<<combinationIndex({5,6,11,12},{0.,0.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[12,70]: "<<combinationIndex({5,6,11,12},{0.,0.,70.,70.},14)<<std::endl;
  std::cout<<"[5,0],[11,0],[12,0],[13,70],[14,70],[15,110]: "<<combinationIndex({5,11,12,13,14,15},{0.,0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[12,0],[13,0],[14,0]: "<<combinationIndex({5,12,13,14},{0.,0.,0.,0.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[8,0],[12,110],[13,110]: "<<combinationIndex({5,6,7,8,12,13},{0.,0.,0.,0.,110.,110.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[8,0],[11,70],[12,110],[13,110]: "<<combinationIndex({6,7,8,11,12,13},{0.,0.,0.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[11,0],[12,0]: "<<combinationIndex({6,7,11,12},{0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[8,0],[11,70],[13,110]: "<<combinationIndex({5,6,7,8,11,13},{0.,0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[12,0],[13,70],[14,70],[15,110]: "<<combinationIndex({5,12,13,14,15},{0.,0.,70.,70.,110.},14)<<std::endl;*/
  //special 5s of new special 6s printout
  /*std::cout<<"[6,0],[11,70],[12,70],[13,110],[14,110]: "<<combinationIndex({6,11,12,13,14},{0.,70.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[11,70],[12,70],[13,110],[14,110]: "<<combinationIndex({5,11,12,13,14},{0.,70.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[12,70],[13,110],[14,110]: "<<combinationIndex({5,6,12,13,14},{0.,0.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[13,110],[14,110]: "<<combinationIndex({5,6,11,13,14},{0.,0.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[12,70],[14,110]: "<<combinationIndex({5,6,11,12,14},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[12,70],[13,110]: "<<combinationIndex({5,6,11,12,13},{0.,0.,70.,70.,110.},14)<<std::endl;

  std::cout<<"[6,0],[11,0],[12,0],[13,0],[14,70]: "<<combinationIndex({6,11,12,13,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[5,0],[11,0],[12,0],[13,0],[14,70]: "<<combinationIndex({5,11,12,13,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[12,0],[13,0],[14,70]: "<<combinationIndex({5,6,12,13,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,0],[13,0],[14,70]: "<<combinationIndex({5,6,11,13,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,0],[12,0],[14,70]: "<<combinationIndex({5,6,11,12,14},{0.,0.,0.,0.,70.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,0],[12,0],[13,0]: "<<combinationIndex({5,6,11,12,13},{0.,0.,0.,0.,0.},14)<<std::endl;

  std::cout<<"[6,0],[8,0],[11,70],[12,110],[13,110]: "<<combinationIndex({6,8,11,12,13},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[8,0],[11,70],[12,110],[13,110]: "<<combinationIndex({5,8,11,12,13},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[11,70],[12,110],[13,110]: "<<combinationIndex({5,6,11,12,13},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[12,110],[13,110]: "<<combinationIndex({5,6,8,12,13},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[11,70],[13,110]: "<<combinationIndex({5,6,8,11,13},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[11,70],[12,110]: "<<combinationIndex({5,6,8,11,12},{0.,0.,0.,70.,70.},14)<<std::endl;

  std::cout<<"[11,0],[12,0],[13,70],[14,70],[15,110]: "<<combinationIndex({11,12,13,14,15},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[12,0],[13,70],[14,70],[15,110]: "<<combinationIndex({5,12,13,14,15},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[11,0],[13,70],[14,70],[15,110]: "<<combinationIndex({5,11,13,14,15},{0.,0.,70.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[11,0],[12,0],[14,70],[15,110]: "<<combinationIndex({5,11,12,14,15},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[11,0],[12,0],[13,70],[15,110]: "<<combinationIndex({5,11,12,13,15},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[11,0],[12,0],[13,70],[14,70]: "<<combinationIndex({5,11,12,13,14},{0.,0.,0.,70.,70.},14)<<std::endl;

  std::cout<<"[6,0],[7,0],[8,0],[12,110],[13,110]: "<<combinationIndex({6,7,8,12,13},{0.,0.,0.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[7,0],[8,0],[12,110],[13,110]: "<<combinationIndex({5,7,8,12,13},{0.,0.,0.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[12,110],[13,110]: "<<combinationIndex({5,6,8,12,13},{0.,0.,0.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[12,110],[13,110]: "<<combinationIndex({5,6,7,12,13},{0.,0.,0.,110.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[8,0],[13,110]: "<<combinationIndex({5,6,7,8,13},{0.,0.,0.,0.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[8,0],[12,110]: "<<combinationIndex({5,6,7,8,12},{0.,0.,0.,0.,110.},14)<<std::endl;

  std::cout<<"[7,0],[8,0],[11,70],[12,110],[13,110]: "<<combinationIndex({7,8,11,12,13},{0.,0.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[6,0],[8,0],[11,70],[12,110],[13,110]: "<<combinationIndex({6,8,11,12,13},{0.,0.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[11,70],[12,110],[13,110]: "<<combinationIndex({6,7,11,12,13},{0.,0.,70.,110.,110.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[8,0],[12,110],[13,110]: "<<combinationIndex({6,7,8,12,13},{0.,0.,0.,110.,110.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[8,0],[11,70],[13,110]: "<<combinationIndex({6,7,8,11,13},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[6,0],[7,0],[8,0],[11,70],[12,110]: "<<combinationIndex({6,7,8,11,12},{0.,0.,0.,70.,110.},14)<<std::endl;

  std::cout<<"[6,0],[7,0],[8,0],[11,70],[13,110]: "<<combinationIndex({6,7,8,11,13},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[7,0],[8,0],[11,70],[13,110]: "<<combinationIndex({5,7,8,11,13},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[8,0],[11,70],[13,110]: "<<combinationIndex({5,6,8,11,13},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[11,70],[13,110]: "<<combinationIndex({5,6,7,11,13},{0.,0.,0.,70.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[8,0],[13,110]: "<<combinationIndex({5,6,7,8,13},{0.,0.,0.,0.,110.},14)<<std::endl;
  std::cout<<"[5,0],[6,0],[7,0],[8,0],[11,70]: "<<combinationIndex({5,6,7,8,11},{0.,0.,0.,0.,70.},14)<<std::endl;*/

  DummyAllowance_=1;
  if(     Index_==2016    || Index_==992     || Index_==1504    || Index_==1760    || Index_==1888    || Index_==1952 || Index_==1984 /*added strange comb*/ || Index_==480 || Index_==864 || Index_==960 || Index_==736 || Index_==928) return 0;
  else if(Index_==3040    || Index_==2528    || Index_==2784    || Index_==2912    || Index_==1976    || Index_==3008 /*added strange comb*/    || Index_==2976 || Index_==2464 || Index_== 2496 || Index_==2400) return 1;
  else if(Index_==6624    || Index_==4576    || Index_==6368    || Index_==6496    || Index_==6560    || Index_==6592 /*added strange comb*/    || Index_==2099424 || Index_==6240    || Index_==12768 || Index_==4320 || Index_==12736 || Index_==12704 || Index_==12640 || Index_==12512 || Index_==8672 || Index_==4576) return 2;
  else if(Index_==2103776 || Index_==2099680 || Index_==2103520 || Index_==2103648 || Index_==2103712 || Index_==2103744 /*added strange comb*/ || Index_==2107872 || Index_==2111936 || Index_==12768 || Index_==2111904 || Index_==10720 || Index_==2103392 || Index_==6297952 || Index_==2107840 || Index_==2107808 || Index_==2107744 || Index_==2107616 || Index_==8672 || Index_==2099680) return 3;
  else if(Index_==2111712 || Index_==2107616 || Index_==12512   || Index_==2111584 || Index_==2111648 || Index_==2111680 /*added strange comb*/ || Index_==6306144 || Index_==2107488 || Index_==2111936 || Index_==6306112 || Index_==6306080 || Index_==6305888 || Index_==4206944 || Index_==2107744 || Index_==2111872 || Index_==2111808 || Index_==2111680 || Index_==12736 || Index_==2107840 || Index_==2103744) return 4;
  else if(Index_==6306016 || Index_==6297824 || Index_==4206816 || Index_==6305888 || Index_==6305952 || Index_==6305984 /*added strange comb*/ || Index_==2111840 || Index_==10464   || Index_==14560 || Index_==14528 || Index_==2272 || Index_==6305856 || Index_==4206688 || Index_==6297696 || Index_==2103488) return 5;
  else if(Index_==4290656 || Index_==4274272 || Index_==4282464 || Index_==92256   || Index_==4223072 || Index_==4290592 || Index_==4290624 /*added strange comb*/ || Index_==4290592 || Index_==6314080|| Index_==6322240 || Index_==6322208 || Index_==4223072 || Index_==2123872 || Index_==6305888) return 6;
  else if(Index_==12679264|| Index_==12662880|| Index_==8480864 || Index_==12611680|| Index_==12679200|| Index_==12679232 /*added strange comb*/|| Index_==6322208 || Index_==14694496 || Index_==28768 || Index_==2127936 || Index_==2119776 || Index_==14710848 || Index_==2127968) return 7;
  else if(Index_==8616032 || Index_==8599648 || Index_==219232  || Index_==8548448 || Index_==8615968 || Index_==8616000) return 8;
  else if(Index_==25393248|| Index_==16996448|| Index_==25258080|| Index_==25325664|| Index_==25393184|| Index_==25393216 /*added strange comb*/|| Index_==6322240 || Index_==6322272 || Index_==2123872 || Index_==473184 || Index_==17199200 || Index_==17266784 || Index_==17266752 || Index_==6314080 || Index_==17131616 || Index_==10512480 || Index_==14710816 || Index_==14710880 || Index_==2127904 || Index_==17266720 || Index_==16996448 || Index_==25393184) return 9;
  else if(Index_==17299488|| Index_==17266720|| Index_==50588   || Index_==17029152|| Index_==17164320|| Index_==17231904 || Index_==17299456 /*added strange comb*/ || Index_==632208 || Index_==25425920 || Index_==6306144 || Index_==25358368|| Index_==25425952|| Index_==25290784 || Index_==17029152 || Index_==8632352) return 10;
  else if(Index_==50853920|| Index_==34060320|| Index_==50583584|| Index_==50718752|| Index_==50853888 /*added strange comb*/|| Index_==25425920 || Index_==25425952 || Index_==25358368 || Index_==8513568 || Index_==251936 || Index_==8581152 || Index_==8648704 || Index_==8648736) return 11;
  else if(Index_==34600992|| Index_==1013792 || Index_==34330656|| Index_==34465824|| Index_==34533408|| Index_==34600960 /*added strange comb*/|| Index_==505888  || Index_==25290784 || Index_==25358368 || Index_==8632352 || Index_==1013760 || Index_==946208) return 12;
  else if(Index_==2095136 || Index_==1554464 || Index_==1824800 || Index_==1959968 || Index_==2027552 || Index_==2095104 /*added strange comb*/ || Index_==50786336 || Index_==9140288 || Index_==34533376 || Index_==2027520 || Index_==1959936 || Index_==1554432 || Index_==1824768) return 13;
  else{
    DummyAllowance_=2; //provide for 4 layer cases
    if(     Index_==480 || Index_==736 || Index_==928 || Index_ ==960) return 0;
    else if(Index_==864) return 1;
    else if(Index_==2272 || Index_==2400 || Index_==2464 || Index_==2496 || Index_==2099648 || Index_==2099616) return 3;
    else if(Index_==4320 || Index_==6240 || Index_==6366 || Index_== 6336 || Index_==2099424 || Index_==2103392 || Index_==2103488 || Index_==6297696 || Index_==2103456 || Index_==6304) return 5;
    else if(Index_==4206688 || Index_==2111520 || Index_==2107488 || Index_==12384 || Index_==2111552) return 7;
    else if(Index_==6305824) return 8;
    else if(Index_==6305856 || Index_==2107488) return 9;
    else if(Index_==1959936 || Index_==1554432 || Index_==2027520 || Index_==1013760 || Index_==1892384 || Index_==1824768 || Index_==1757216 || Index_==946208 || Index_==1486880) return 13;
  }
  return 999;
}

IndexInfo PDDS::IndexDefiner(unsigned int target_){ //redefines the compact index info to the target combination kind
  IndexInfo in_;
  if(target_==0){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=7;
    in_.layers[3]=8;
    in_.layers[4]=9;
    in_.layers[5]=10;
  }
  else if(target_==1){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=7;
    in_.layers[3]=8;
    in_.layers[4]=9;
    in_.layers[5]=11;
  }
  else if(target_==2){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=7;
    in_.layers[3]=8;
    in_.layers[4]=11;
    in_.layers[5]=12;
  }
  else if(target_==3){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=7;
    in_.layers[3]=8;
    in_.layers[4]=11;
    in_.layers[5]=12;
  }
  else if(target_==4){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=7;
    in_.layers[3]=11;
    in_.layers[4]=12;
    in_.layers[5]=13;
  }
  else if(target_==5){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=7;
    in_.layers[3]=11;
    in_.layers[4]=12;
    in_.layers[5]=13;
  }
  else if(target_==6){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=11;
    in_.layers[3]=12;
    in_.layers[4]=13;
    in_.layers[5]=14;
  }
  else if(target_==7){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=11;
    in_.layers[3]=12;
    in_.layers[4]=13;
    in_.layers[5]=14;
  }
  else if(target_==8){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=11;
    in_.layers[3]=12;
    in_.layers[4]=13;
    in_.layers[5]=14;
  }
  else if(target_==9){
    in_.layers[0]=5;
    in_.layers[1]=6;
    in_.layers[2]=11;
    in_.layers[3]=12;
    in_.layers[4]=13;
    in_.layers[5]=14;
  }
  else if(target_==10){
    in_.layers[0]=5;
    in_.layers[1]=11;
    in_.layers[2]=12;
    in_.layers[3]=13;
    in_.layers[4]=14;
    in_.layers[5]=15;
  }
  else if(target_==11){
    in_.layers[0]=5;
    in_.layers[1]=11;
    in_.layers[2]=12;
    in_.layers[3]=13;
    in_.layers[4]=14;
    in_.layers[5]=15;
  }
  else if(target_==12){
    in_.layers[0]=5;
    in_.layers[1]=11;
    in_.layers[2]=12;
    in_.layers[3]=13;
    in_.layers[4]=14;
    in_.layers[5]=15;
  }
  else if(target_==13){
    in_.layers[0]=5;
    in_.layers[1]=11;
    in_.layers[2]=12;
    in_.layers[3]=13;
    in_.layers[4]=14;
    in_.layers[5]=15;
  }
  return in_;
}

unsigned int PDDS::CombinationKind(std::vector<unsigned> layer_, std::vector<float> radii_){ //gives back the cut region depending on the layers and radii of incoming stubs
  std::vector<int> Layers;
  std::vector<double> Radii;
  bool set[11]={false,false,false,false,false,false,false,false,false,false,false};
  for(unsigned s=0; s<layer_.size(); ++s){
    const unsigned layer=layer_[s];
    if(layer==5 && !set[0])                      {Layers.push_back(layer); Radii.push_back(radii_[s]); set[0]=true;}
    else if(layer==6 && !set[1])		 {Layers.push_back(layer); Radii.push_back(radii_[s]); set[1]=true;}
    else if(layer==7 && !set[2])		 {Layers.push_back(layer); Radii.push_back(radii_[s]); set[2]=true;}
    else if(layer==8 && !set[3])		 {Layers.push_back(layer); Radii.push_back(radii_[s]); set[3]=true;}
    else if(layer==9 && !set[4])		 {Layers.push_back(layer); Radii.push_back(radii_[s]); set[4]=true;}
    else if(layer==10 && !set[5])		 {Layers.push_back(layer); Radii.push_back(radii_[s]); set[5]=true;}
    else if((layer==11 || layer==18) && !set[6]) {Layers.push_back(11); Radii.push_back(radii_[s]); set[6]=true;}
    else if((layer==12 || layer==19) && !set[7]) {Layers.push_back(12); Radii.push_back(radii_[s]); set[7]=true;}
    else if((layer==13 || layer==20) && !set[8]) {Layers.push_back(13); Radii.push_back(radii_[s]); set[8]=true;}
    else if((layer==14 || layer==21) && !set[9]) {Layers.push_back(14); Radii.push_back(radii_[s]); set[9]=true;}
    else if((layer==15 || layer==22) && !set[10]){Layers.push_back(15); Radii.push_back(radii_[s]); set[10]=true;}
  }
  const unsigned int output=IndexInterpreter(combinationIndex(Layers,Radii,14));
  if(verbose_){
    std::cout<<"MarcoIndex "<<combinationIndex(Layers,Radii,14)<<"-> category "<<output<<std::endl;
    for(unsigned i=0; i<layer_.size(); ++i) std::cout<<"["<<layer_[i]<<","<<radii_[i]<<"], ";
    std::cout<<std::endl;
  }
  /*if(output==999){
    std::cout<<combinationIndex(Layers,Radii,14)<<"->"<<output<<std::endl;
    for(unsigned i=0; i<layer_.size(); ++i) std::cout<<"["<<layer_[i]<<","<<radii_[i]<<"], ";
    std::cout<<std::endl;
  }*/
  return output;
}

void PDDS::initialize(std::vector<unsigned> layer, std::vector<float> radii){
  Index=PDDS::CombinationKind(layer,radii);
  cuts_[0]=99.; cuts_[1]=99.; cuts_[2]=99.; //no cuts on 
  if(Index==0)      {cuts_[0]=3.0;cuts_[1]=3.0;cuts_[2]=2.0;}
  else if(Index==1) {cuts_[0]=2.0;cuts_[1]=2.5;cuts_[2]=2.0;}
  else if(Index==2) {cuts_[0]=2.0;cuts_[1]=2.5;cuts_[2]=2.0;}
  else if(Index==3) {cuts_[0]=2.0;cuts_[1]=3.0;cuts_[2]=2.0;}
  else if(Index==4) {cuts_[0]=2.0;cuts_[1]=3.0;cuts_[2]=2.0;}
  else if(Index==5) {cuts_[0]=3.0;cuts_[1]=3.5;cuts_[2]=2.0;}
  else if(Index==6) {cuts_[0]=2.0;cuts_[1]=2.0;cuts_[2]=2.0;}
  else if(Index==7) {cuts_[0]=3.0;cuts_[1]=3.5;cuts_[2]=3.0;}
  else if(Index==8) {cuts_[0]=3.0;cuts_[1]=2.5;cuts_[2]=3.5;}
  else if(Index==9) {cuts_[0]=2.0;cuts_[1]=2.0;cuts_[2]=2.0;}
  else if(Index==10){cuts_[0]=3.0;cuts_[1]=2.5;cuts_[2]=3.5;}
  else if(Index==11){cuts_[0]=3.0;cuts_[1]=2.0;cuts_[2]=3.5;}
  else if(Index==12){cuts_[0]=2.0;cuts_[1]=3.0;cuts_[2]=2.0;}
  else if(Index==13){cuts_[0]=2.0;cuts_[1]=2.0;cuts_[2]=2.0;}
}

std::vector<PairAssignment> PDDS::Stage1Pair(std::vector<std::pair<unsigned, float> > first, std::vector<std::pair<unsigned, float> > second, float cut) //applies pairwise PDDS cuts and organizes layer-pairs
{  
  //initialize pair assignment struct, loop and filter pairs
  std::vector<PairAssignment> PairedOnes;
  for(unsigned i=0; i<first.size(); ++i){
    for(unsigned j=0; j<second.size(); ++j){
      bool Survivor=true;
      //if(first[i].first==999999999 && second[j].first==999999999) Survivor=false; //catch double dummy cases
      /*else */if(fabs(first[i].second-second[j].second)>cut && !(first[i].first==999999999 || second[j].first==999999999) ) Survivor=false; //pairwise DDS cut, spares any pair with only a single dummy
      PairedOnes.push_back(PairAssignment(first[i].first,second[j].first,Survivor));
    }
  }
  //give back assignment map
  return PairedOnes;
}

//combination(stubAddresses,pass)
std::vector<std::vector<unsigned> > PDDS::CombinationBuilder(std::vector<PairAssignment> Layer10, std::vector<PairAssignment> Layer32, std::vector<PairAssignment> Layer54)
{
  //initialize return object
  std::vector<std::vector<unsigned> >  Combinations;
  if(verbose_) std::cout<<"----- Running PDDS with cuts ["<<cuts_[0]<<","<<cuts_[1]<<","<<cuts_[2]<<"] ------"<<std::endl;
  //loop over all permutations
  for(unsigned i=0; i<Layer10.size(); ++i){
    for(unsigned j=0; j<Layer32.size(); ++j){
      for(unsigned k=0; k<Layer54.size(); ++k){
	std::vector<unsigned> layers={Layer10[i].Stub2Address, Layer10[i].Stub1Address, Layer32[j].Stub2Address, Layer32[j].Stub1Address, Layer54[k].Stub2Address, Layer54[k].Stub1Address}; //combination

	//catch >1 dummy cases, always put dummy in the first position per layer!
	unsigned dummyCounter=0;
	for(unsigned l=0; l<layers.size(); ++l) if(layers[l]==999999999) ++dummyCounter;
	if(dummyCounter>DummyAllowance_){
	  //if(verbose_) std::cout<<"Failed DC "<<Layer10[i].Stub2Address<<", "<<Layer10[i].Stub1Address<<", "<<Layer32[j].Stub2Address<<", "<<Layer32[j].Stub1Address<<", "<<Layer54[k].Stub2Address<<", "<<Layer54[k].Stub1Address<<std::endl;
	  continue;
	}

	//check for DDS failure
	if(!(Layer10[i].SurvivingPair * Layer32[j].SurvivingPair * Layer54[k].SurvivingPair)){
	  if(verbose_) std::cout<<"Failed DS: "<<Layer10[i].Stub2Address<<", "<<Layer10[i].Stub1Address<<", "<<Layer32[j].Stub2Address<<", "<<Layer32[j].Stub1Address<<", "<<Layer54[k].Stub2Address<<", "<<Layer54[k].Stub1Address<<std::endl;
	  continue; //3-logic part
	}
	
	if(verbose_) std::cout<<"Passed "<<Layer10[i].Stub2Address<<", "<<Layer10[i].Stub1Address<<", "<<Layer32[j].Stub2Address<<", "<<Layer32[j].Stub1Address<<", "<<Layer54[k].Stub2Address<<", "<<Layer54[k].Stub1Address<<std::endl;
	//build combination object
	Combinations.push_back(layers);
      }
    }
  }
  return Combinations;
}

float PDDS::SebCoarsener(float DS_, int layer_, float z_, unsigned ring_){
  if(layer_==5 || layer_==6 || layer_==7){
    if(fabs(DS_)<=0.5)        return 0.;
    else if (fabs(DS_)==1.)   return sgn(DS_)*1.;
    else if (fabs(DS_)<=2.)   return sgn(DS_)*2.;
    else return sgn(DS_)*3.;
  }
  else if(layer_==8){
    if(fabs(DS_)<=0.5)        return 0.;
    else if(fabs(DS_)==1.)    return sgn(DS_)*1.;
    else if(fabs(DS_)<=2.)    return sgn(DS_)*2.;
    else return DS_;
  }
  else if(layer_==9){
    if(fabs(DS_)<=0.5)        return 0.;
    else if(fabs(DS_)==1.)    return sgn(DS_)*1.;
    else if(fabs(DS_)<=2.5)   return sgn(DS_)*2.;
    else if(fabs(DS_)<=3.5)   return sgn(DS_)*3.;
    else return DS_;
  }
  else if(layer_==10){
    if(fabs(DS_)<=0.5)        return 0.;
    else if(fabs(DS_)==1.)    return sgn(DS_)*1.;
    else if(fabs(DS_)<=2.5)   return sgn(DS_)*2.;
    else if(fabs(DS_)<=4)     return sgn(DS_)*3.;
    else if(fabs(DS_)<=5)     return sgn(DS_)*4.;
    else return DS_;
  }
  else if(layer_>=11){
    if(ring_<4){
      if(fabs(DS_)<=0.5)      return 0.;
      else if(fabs(DS_)==1.)  return -1.*sgn(z_)*sgn(DS_)*1.;
      else if(fabs(DS_)==1.5) return -1.*sgn(z_)*sgn(DS_)*2.;
      else if(fabs(DS_)==2.)  return -1.*sgn(z_)*sgn(DS_)*3.;
    }
    else if(ring_<7){
      if(fabs(DS_)<=1.)       return 0.;
      else if(fabs(DS_)==1.5) return -1.*sgn(z_)*sgn(DS_)*1.;
      else if(fabs(DS_)==2.)  return -1.*sgn(z_)*sgn(DS_)*2.;
      else if(fabs(DS_)==2.5) return -1.*sgn(z_)*sgn(DS_)*3.;
    }
    else if(ring_==7){
      if(fabs(DS_)<=1.)       return 0.;
      else if(fabs(DS_)<=2.)  return -1.*sgn(z_)*sgn(DS_)*1.;
      else if(fabs(DS_)==2.5) return -1.*sgn(z_)*sgn(DS_)*2.;
      else if(fabs(DS_)==3.)  return -1.*sgn(z_)*sgn(DS_)*3.;
    }
    else if(ring_==8){
      if(fabs(DS_)<=1.)       return 0.;
      else if(fabs(DS_)<=2.)  return -1.*sgn(z_)*sgn(DS_)*1.;
      else if(fabs(DS_)<=3.)  return -1.*sgn(z_)*sgn(DS_)*2.;
      else if(fabs(DS_)==3.5) return -1.*sgn(z_)*sgn(DS_)*3.;
    }
    else{
      if(fabs(DS_)<=0.5)      return 0.;
      else if(fabs(DS_)==1.)  return -1.*sgn(z_)*sgn(DS_)*1.;
      else if(fabs(DS_)<=2.)  return -1.*sgn(z_)*sgn(DS_)*2.;
      else if(fabs(DS_)<=3.)  return -1.*sgn(z_)*sgn(DS_)*3.;     
      else if(fabs(DS_)<=4.)  return -1.*sgn(z_)*sgn(DS_)*4.;   
      else if(fabs(DS_)==4.5) return -1.*sgn(z_)*sgn(DS_)*5.;   
      else if(fabs(DS_)==5.)  return -1.*sgn(z_)*sgn(DS_)*6.;
      else if(fabs(DS_)==5.5) return -1.*sgn(z_)*sgn(DS_)*7.; 
    }
  }
  return -999.;
}

//get the basic combination layout
std::vector<std::vector<unsigned> > PDDS::CombinationAll(std::vector<Stub> stubs_){
  std::vector<std::vector<unsigned> > CombinationBase_={{},{},{},{},{},{},{},{},{},{},{}};
  for(unsigned s=0; s<stubs_.size(); ++s){
    const unsigned stub=s;
    const unsigned layer=stubs_[s].layer();
    if(layer==5) CombinationBase_[0].push_back(stub);
    else if(layer==6)  CombinationBase_[1].push_back(stub);
    else if(layer==7)  CombinationBase_[2].push_back(stub);
    else if(layer==8)  CombinationBase_[3].push_back(stub);
    else if(layer==9)  CombinationBase_[4].push_back(stub);
    else if(layer==10) CombinationBase_[5].push_back(stub);
    else if(layer==11 || layer==18) CombinationBase_[6].push_back(stub);
    else if(layer==12 || layer==19) CombinationBase_[7].push_back(stub);
    else if(layer==13 || layer==20) CombinationBase_[8].push_back(stub);
    else if(layer==14 || layer==21) CombinationBase_[9].push_back(stub);
    else if(layer==15 || layer==22) CombinationBase_[10].push_back(stub);
  }
  for(unsigned i=0; i<CombinationBase_.size(); ++i){//at most 4 stubs per layer cutoff
    if(CombinationBase_[i].size()>4){
      std::vector<unsigned> replacement;
      for(unsigned j=CombinationBase_[i].size()-4; j<CombinationBase_[i].size(); ++j) replacement.push_back(CombinationBase_[i][j]);
      CombinationBase_[i]=replacement;
      if(CombinationBase_[i].size()!=4) std::cout<<"didn't work! Still has "<<CombinationBase_[i].size()<<" entries!!!"<<std::endl;
    }
  }
  return CombinationBase_;
}

std::vector<StubsCombination> PDDS::combine(Road & road){

  //option to activate detailed printout
  verbose_=false;

  std::vector<StubsCombination> Result={};

  //find out the road index
  std::vector<unsigned> layer;
  std::vector<float> radii;
  std::vector<std::vector<Stub> > Stubs = road.getStubs();
  std::vector<Stub> FlatStubs;
  if(verbose_) std::cout<<"_____incoming new road_____"<<std::endl;
  for(unsigned l=0; l<Stubs.size(); ++l){
    if(verbose_) std::cout<<"Layer "<<l<<": ";
    if(Stubs[l].size()>0){
      layer.push_back(Stubs[l][0].layer());
      radii.push_back(Stubs[l][0].R());
    }
    for(unsigned s=0; s<Stubs[l].size(); ++s){
      FlatStubs.push_back(Stubs[l][s]);
      if(verbose_) std::cout<<"["<<FlatStubs.size()-1<<","<<Stubs[l][s].DeltaS()<<"],";
    }
    if(verbose_) std::cout<<std::endl;
  }

  PDDS::initialize(layer,radii);
  IndexInfo info = PDDS::IndexDefiner(Index);
  if(Index==999 && verbose_){//exception for unknown combination!!!
    std::cout<<"PDDS failed!!! Do not know combination with index "<<Index<<std::endl;
    return Result;
  }
  std::vector<std::vector<unsigned> > CombinationAll_=PDDS::CombinationAll(FlatStubs);

  if(SetFive_){//insert dummies in every layer
    CombinationAll_[info.layers[0]-5].insert(CombinationAll_[info.layers[0]-5].begin(),999999999);
    CombinationAll_[info.layers[1]-5].insert(CombinationAll_[info.layers[1]-5].begin(),999999999);
    CombinationAll_[info.layers[2]-5].insert(CombinationAll_[info.layers[2]-5].begin(),999999999);
    CombinationAll_[info.layers[3]-5].insert(CombinationAll_[info.layers[3]-5].begin(),999999999);
    CombinationAll_[info.layers[4]-5].insert(CombinationAll_[info.layers[4]-5].begin(),999999999);
    CombinationAll_[info.layers[5]-5].insert(CombinationAll_[info.layers[5]-5].begin(),999999999);
  }
  else{//insert dummies for empty layer
    if(!CombinationAll_[info.layers[0]-5].size()) CombinationAll_[info.layers[0]-5].insert(CombinationAll_[info.layers[0]-5].begin(),999999999);
    if(!CombinationAll_[info.layers[1]-5].size()) CombinationAll_[info.layers[1]-5].insert(CombinationAll_[info.layers[1]-5].begin(),999999999);
    if(!CombinationAll_[info.layers[2]-5].size()) CombinationAll_[info.layers[2]-5].insert(CombinationAll_[info.layers[2]-5].begin(),999999999);
    if(!CombinationAll_[info.layers[3]-5].size()) CombinationAll_[info.layers[3]-5].insert(CombinationAll_[info.layers[3]-5].begin(),999999999);
    if(!CombinationAll_[info.layers[4]-5].size()) CombinationAll_[info.layers[4]-5].insert(CombinationAll_[info.layers[4]-5].begin(),999999999);
    if(!CombinationAll_[info.layers[5]-5].size()) CombinationAll_[info.layers[5]-5].insert(CombinationAll_[info.layers[5]-5].begin(),999999999);
  }
  
  //fill layer-ordered stub containers for DDS processing
  std::vector<std::vector<std::pair<unsigned,float> > > Layer;
  for(unsigned l=0; l<6; ++l){
    const unsigned layer=info.layers[l];
    std::vector<std::pair<unsigned,float> > OneLayer;
    Layer.push_back(OneLayer);
    for(unsigned s=0; s<CombinationAll_[layer-5].size(); ++s){
      const unsigned stub=CombinationAll_[layer-5][s];
      float DSinit=0.;
      if(stub!=999999999){
	DSinit=FlatStubs[stub].DeltaS();
	Layer[l].push_back(std::make_pair(stub,PDDS::SebCoarsener(DSinit,layer,FlatStubs[stub].z(),FlatStubs[stub].ring()))); //use coarsening for DDS
      }
      else{
	Layer[l].push_back(std::make_pair(stub,DSinit));
      }
    }
  }
  //build reduced layer pairs
  std::vector<PairAssignment> Layer01=PDDS::Stage1Pair(Layer[1],Layer[0],cuts_[0]);
  std::vector<PairAssignment> Layer23=PDDS::Stage1Pair(Layer[3],Layer[2],cuts_[1]);
  std::vector<PairAssignment> Layer45=PDDS::Stage1Pair(Layer[5],Layer[4],cuts_[2]);
  //contain 
  std::vector<std::vector<unsigned> > combinations=PDDS::CombinationBuilder(Layer01,Layer23,Layer45);
  //convert combinations into StubCombinations
  if(verbose_) std::cout<<"-_-_-_-_-_Pushing out final combinations_-_-_-_-_-"<<std::endl;
  for(unsigned c=0; c<combinations.size(); ++c){
    StubsCombination element;
    if(verbose_) std::cout<<"combination "<<c<<": ";
    for(unsigned s=0; s<combinations[c].size(); ++s){
      const unsigned address = combinations[c][s];
      if(address!=999999999){
	element.pushStub(FlatStubs[address]);
	if(verbose_) std::cout<<address<<", ";
      }
    }
    if(verbose_) std::cout<<std::endl;
    //fillGenInfo(element); //does not do anything
    Result.push_back(element);
  }
  return Result;
}
