#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include<utility>      // std::pair
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSPUenergy.hh"
#include "HGCSSSimpleHit.hh"

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

double DeltaR(double eta1,double phi1,double eta2,double phi2){
  double dr=99999.;
  double deta=fabs(eta1-eta2);
  double dphi=fabs(phi1-phi2);
  if(dphi>TMath::Pi()) dphi=2.*TMath::Pi()-dphi;
  dr=sqrt(deta*deta+dphi*dphi);
  return dr;
}


int main(int argc, char** argv){//main  


  //Input output and config options
  std::string cfg;
  unsigned pNevts;
  //std::string inFilePath;
  std::string outFilePath;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
  unsigned debug;
  //double etamean;
  //double deta;
  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("digifilePath", po::value<std::string>(&digifilePath)->default_value(""))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    //("etamean,e",      po::value<double>(&etamean)->default_value(2.8))
    //("deta",      po::value<double>(&deta)->default_value(0.05))
    ;
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);


  std::string inFilePath = filePath+simFileName;

  std::cout << " -- Input parameters: " << std::endl
            << " -- Input file path: " << filePath << std::endl
            << " -- Digi Input file path: " << digifilePath << std::endl
    	    << " -- Output file path: " << outFilePath << std::endl
    //	    << " -- mean eta: " << etamean 
	    << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events per run." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;


  //
  // hardcoded
  //


  //global threshold to reduce size of noise hits
  const double threshMin = 0.5;

  std::cout << " ---- Selection settings: ---- " << std::endl
	    << " -------threshMin " << threshMin << std::endl
	    << " ------------------------------------------" << std::endl;



  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////


  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
  else
    inputrec << digifilePath << "/" << recoFileName;

  std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;

  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos)
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
    else {
      std::cout << " -- Error in getting information from simfile!" << std::endl;
      return 1;
    }
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstrsim;
      std::ostringstream lstrrec;
      lstrsim << inputsim.str() << "_run" << i << ".root";
      if (testInputFile(lstrsim.str(),simFile)){
        if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
        else {
	  std::cout << " -- Error in getting information from simfile!" << std::endl;
          return 1;
        }
      }
      else continue;
      lstrrec << inputrec.str() << "_run" << i << ".root";
      if (!testInputFile(lstrrec.str(),recFile)) continue;
      lSimTree->AddFile(lstrsim.str().c_str());
      lRecTree->AddFile(lstrrec.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  //assert(info);

  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  const unsigned shape = info->shape();
  const double cellSize = info->cellSize();
  const double calorSizeXY = info->calorSizeXY();

  bool isTBsetup = (model != 2);
  bool bypassR = false;
  if (isTBsetup) bypassR = true;

  HGCSSDetector & myDetector = theDetector();
  myDetector.buildDetector(versionNumber,true,false,bypassR);

  //corrected for Si-Scint overlap
  const unsigned nLayers = 52;//


  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model << std::endl
	    << " -- cellSize = " << cellSize
	    << ", shape = " << shape
	    << ", nLayers = " << nLayers
	    << std::endl;
  HGCSSGeometryConversion geomConv(model,cellSize,bypassR,3);

  geomConv.setXYwidth(calorSizeXY);
  geomConv.setVersion(versionNumber);
  
  if (shape==2) geomConv.initialiseDiamondMap(calorSizeXY,10.);
  else if (shape==3) geomConv.initialiseTriangleMap(calorSizeXY,10.*sqrt(2.));
  else if (shape==1) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  else if (shape==4) geomConv.initialiseSquareMap(calorSizeXY,10.);

  //square map for BHCAL
  //geomConv.initialiseSquareMap1(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),0.01745);//eta phi segmentation
  //geomConv.initialiseSquareMap2(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),0.02182);//eta phi segmentation
  geomConv.initialiseSquareMap1(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),TMath::Pi()*2./360.);//eta phi segmentation
  geomConv.initialiseSquareMap2(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),TMath::Pi()*2./288.);//eta phi segmentation
  std::vector<unsigned> granularity;
  granularity.resize(myDetector.nLayers(),1);
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();




  std::ostringstream label;
  //root doesn't like . in branch names.....
  //label << "hits";
  //miptree->Branch(label.str().c_str(),"std::vector<HGCSSSimpleHit>",&miphitvec);

  double energy_max=-1.;

  TH1F* h_energy = new TH1F("h_energy","hit energy",1000,0.,5.);
  TH1F* h_z = new TH1F("h_z","z of hit",5000,3100.,5200);
  TH1F* h_z1 = new TH1F("h_z1","z of hit",5000,3150.,3550);
  TH1F* h_z2 = new TH1F("h_z2","z of hit",5000,3550.,5200);
  TH2F* h_xy = new TH2F("h_xy","xy of hit",1000,-2000,2000,1000,-2000,2000);
  TH2F* h_etaphi = new TH2F("h_etaphi","etaphi of hit",1000,1,3.5,1000,-7,7);
  TH2F* h_getaphi = new TH2F("h_getaphi","gen part etaphi of hit",1000,1,3.5,1000,-7,7);
  TH1F* h_l = new TH1F("h_l","layer of hit",80,0.,80.);
  TH1F* h_l2 = new TH1F("h_l2","layer of hit",30,50,80.);
  TH2F* h_zl = new TH2F("h_zl","z vs l of hit",5000,4300.,5200,25,30.,55.);
  /////////////////////////
  //Bryans analysis stuff//
  /////////////////////////
  TProfile* h_el = new TProfile("h_el","energy per layer",80,0,80,0,100);
 
  TH2Poly* map_1 = new TH2Poly();
  //map_1->Honeycomb(-2803.17,-2790.5,6.49635,575,497);
  //map_1 = geomConv.hexagonMap();
  map_1 = geomConv.hexagonMap();
  map_1->SetTitle("map_1");
  map_1->SetName("map_1");

  //TH2Poly* map_1_0 = (TH2Poly*) map_1->Clone("map_1_0"); map_1_0->SetTitle("map_1_0");  map_1_0->SetName("map_1_0"); 
  //TH2Poly* map_1_1 = (TH2Poly*) map_1->Clone("map_1_1"); map_1_1->SetTitle("map_1_1");  map_1_1->SetName("map_1_1");

  char title[100];
  TH2Poly* map_1_layer[52];
  bool develop=false;
  if (develop){
  for (int ilayer=36;ilayer<=51;ilayer++){
    sprintf(title,"map_1_%d",ilayer);
    //map_1_layer[ilayer] = new TH2Poly();
    //map_1_layer[ilayer] = geomConv.hexagonMap();
    map_1_layer[ilayer] = (TH2Poly*)map_1->Clone(title);
    map_1_layer[ilayer]->SetTitle(title);
    map_1_layer[ilayer]->SetName(title);
  }
  }

  //TH2Poly* map_2 = new TH2Poly();
  //map_2 = geomConv.squareMap1();
  //map_2 = geomConv.squareMap1();

  //TH2Poly* map_3 = new TH2Poly();
  //map_3 = geomConv.squareMap2();
  //map_3 = geomConv.squareMap2();

  //---- xmin = 1.4, ymin = -3.14159 side = 0.01745, nx = 91, ny=360
  //---- xmin = 1.4, ymin = -3.14159 side = 0.02182, nx = 73, ny=287 --> 288
  double rbins2[92];
  double rbins3[74];
  TH2F* map_TH2F_2[4];
  TH2F* map_TH2F_3[12];
  double z_layer[69]={
  3198.0,    3207.1,    3222.4,    3231.5,    3246.8,
  3255.9,    3271.2,    3280.3,    3295.6,    3304.7,
  3320.0,    3329.1,    3344.4,    3353.5,    3368.8,
  3377.9,    3393.2,    3402.3,    3417.6,    3426.7,
  3442.0,    3451.1,    3466.4,    3475.5,    3490.8,
  3499.9,    3515.2,    3524.3,    3577.4,    3626.4,
  3675.4,    3724.4,    3773.4,    3822.4,    3871.4,
  3920.4,    3969.4,    4020.3,    4071.2,    4122.1,
  4206.0,    4289.9,    4373.8,    4457.7,    4541.6,
  4625.5,    4709.4,    4793.3,    4877.2,    4961.1,
  5045.0,    5128.9,       0.0,    3971.2,    4022.1,
  4073.0,    4123.9,    4207.8,    4291.7,    4375.6,
  4459.5,    4543.4,    4627.3,    4711.2,    4795.1,
  4879.0,    4962.9,    5046.8,    5130.7};

  int eta_index[69]={
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,0,0, 0,0,0,0,0,
    0,0,0,19,24, 25,30,29,32,37,
    41,41,42,43,44, 45,45,46,47};

  int bin_exclude=2;
  double r_cut[69]; // r_cut values to exclude a few boundary hits

  // BH fine part
  for (int ilayer=36;ilayer<=39;ilayer++){
    int ilayer_org=ilayer+17;
    for (int ibin=0;ibin<=91;ibin++){
      double eta=1.4+double(ibin)*0.01745;
      double z = z_layer[ilayer_org];
      rbins2[91-ibin]=z*tan(2.*atan(exp(-eta)));    
    }    
    sprintf(title,"map_TH2F_2_%d",ilayer);
    map_TH2F_2[ilayer-36] = new TH2F(title,title,360,-1.*TMath::Pi(),TMath::Pi(),91,rbins2);
    double eta_tmp=1.4+double(eta_index[ilayer_org]-bin_exclude)*0.01745;
    double z_tmp = z_layer[ilayer_org];
    r_cut[ilayer_org] = z_tmp*tan(2.*atan(exp(-eta_tmp)));    
  }  
  // BH coarse part
  for (int ilayer=40;ilayer<=51;ilayer++){
    int ilayer_org=ilayer+17;
    for (int ibin=0;ibin<=73;ibin++){
      double eta=1.4+double(ibin)*0.02182;
      double z = z_layer[ilayer_org];
      rbins3[73-ibin]=z*tan(2.*atan(exp(-eta)));    
    }    
    sprintf(title,"map_TH2F_3_%d",ilayer);
    map_TH2F_3[ilayer-40] = new TH2F(title,title,288,-1.*TMath::Pi(),TMath::Pi(),73,rbins3);
    double eta_tmp=1.4+double(eta_index[ilayer_org]-bin_exclude)*0.02182;
    double z_tmp = z_layer[ilayer_org];
    r_cut[ilayer_org] = z_tmp*tan(2.*atan(exp(-eta_tmp)));    
  } 

  /////////////////////////
  TH2F* h_zx = new TH2F("h_zx","zx of hit",5000,3100,5200,1000,-2000,2000);
  TH2F* h_zx10000 = new TH2F("h_zx10000","zx of hit",10000,3100,5200,1000,-2000,2000);
  TH2F* h_zx1000 = new TH2F("h_zx1000","zx of hit",1000,3100,5200,1000,-2000,2000);
  TH3F* h_xyz = new TH3F("h_xyz","xyz of hit",100,-2000.,2000.,100,-2000.,2000.,500,3100.,5200.); // 3d histo
  TH2F* h_yz = new TH2F("h_yz","yz of hit",5000,3100,5200,1000,-2000,2000);
  TH2F* h_zx1 = new TH2F("h_zx1","zx of hit",5000,3150,3550,4000,-2000,2000);
  TH2F* h_yz1 = new TH2F("h_yz1","yz of hit",5000,3150,3550,4000,-2000,2000);
  TH2F* h_zx2 = new TH2F("h_zx2","zx of hit",5000,3550,5200,4000,-2000,2000);
  TH2F* h_yz2 = new TH2F("h_yz2","yz of hit",5000,3550,5200,4000,-2000,2000);
  TH2F* h_zr = new TH2F("h_zr","zr of hit",5000,3100,5200,1000,0,2000);
  
  //////layer histos///////

  TH2F* h_nszxl = new TH2F("h_nszxl","zx of hit not scint (layers)",1400,3800,5200,3000,-2000,2000);
  TH2F* h_szxl = new TH2F("h_szxl","zx of hit scint (layers)",1400,3800,5200,3000,-2000,2000);
  TH2F* h_nszx = new TH2F("h_nszx","zx of hit not scint",1900,3100,5200,4000,-2000,2000);
  TH2F* h_szx = new TH2F("h_szx","zx of hit scint",1900,3100,5200,4000,-2000,2000);
  TH2F* h_nsxy = new TH2F("h_nsxy","xy of hit not scint",2000,-2000,2000,2000,-2000,2000);
  TH2F* h_sxy = new TH2F("h_sxy","xy of hit scint",2000,-2000,2000,2000,-2000,2000);
  TH2F* h_nsxyl = new TH2F("h_nsxyl","xy of hit not scint (layers)",3000,-2000,2000,3000,-2000,2000);
  TH2F* h_sxyl = new TH2F("h_sxyl","xy of hit scint (layers)",3000,-2000,2000,3000,-2000,2000);
    
  TH2F* h_nszx36 = new TH2F("h_nszx36","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx37 = new TH2F("h_nszx37","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx38 = new TH2F("h_nszx38","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx39 = new TH2F("h_nszx39","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx40 = new TH2F("h_nszx40","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx41 = new TH2F("h_nszx41","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx42 = new TH2F("h_nszx42","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx43 = new TH2F("h_nszx43","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx44 = new TH2F("h_nszx44","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx45 = new TH2F("h_nszx45","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx46 = new TH2F("h_nszx46","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx47 = new TH2F("h_nszx47","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx48 = new TH2F("h_nszx48","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx49 = new TH2F("h_nszx49","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx50 = new TH2F("h_nszx50","zx of hit not scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_nszx51 = new TH2F("h_nszx51","zx of hit not scint",5000,3100,5200,1000,-1200,1200);

  TH2F* h_szx36 = new TH2F("h_szx36","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx37 = new TH2F("h_szx37","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx38 = new TH2F("h_szx38","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx39 = new TH2F("h_szx39","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx40 = new TH2F("h_szx40","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx41 = new TH2F("h_szx41","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx42 = new TH2F("h_szx42","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx43 = new TH2F("h_szx43","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx44 = new TH2F("h_szx44","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx45 = new TH2F("h_szx45","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx46 = new TH2F("h_szx46","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx47 = new TH2F("h_szx47","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx48 = new TH2F("h_szx48","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx49 = new TH2F("h_szx49","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx50 = new TH2F("h_szx50","zx of hit scint",5000,3100,5200,1000,-1200,1200);
  TH2F* h_szx51 = new TH2F("h_szx51","zx of hit scint",5000,3100,5200,1000,-1200,1200);

  
  ///////////END///////////

  TH2F* h_nsxy36 = new TH2F("h_nsxy36","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy37 = new TH2F("h_nsxy37","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy38 = new TH2F("h_nsxy38","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy39 = new TH2F("h_nsxy39","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy40 = new TH2F("h_nsxy40","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy41 = new TH2F("h_nsxy41","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy42 = new TH2F("h_nsxy42","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy43 = new TH2F("h_nsxy43","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy44 = new TH2F("h_nsxy44","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy45 = new TH2F("h_nsxy45","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy46 = new TH2F("h_nsxy46","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy47 = new TH2F("h_nsxy47","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy48 = new TH2F("h_nsxy48","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy49 = new TH2F("h_nsxy49","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy50 = new TH2F("h_nsxy50","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy51 = new TH2F("h_nsxy51","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);

  TH2F* h_sxy36 = new TH2F("h_sxy36","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy37 = new TH2F("h_sxy37","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy38 = new TH2F("h_sxy38","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy39 = new TH2F("h_sxy39","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy40 = new TH2F("h_sxy40","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy41 = new TH2F("h_sxy41","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy42 = new TH2F("h_sxy42","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy43 = new TH2F("h_sxy43","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy44 = new TH2F("h_sxy44","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy45 = new TH2F("h_sxy45","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy46 = new TH2F("h_sxy46","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy47 = new TH2F("h_sxy47","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy48 = new TH2F("h_sxy48","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy49 = new TH2F("h_sxy49","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy50 = new TH2F("h_sxy50","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy51 = new TH2F("h_sxy51","xy of hit scint",1000,-1200,1200,1000,-1200,1200);

  TH2F* h_Egenreco = new TH2F("h_Egenreco","E reco sum versus gen",1000,0.,1000.,100,0.,20.);
  TH1F* h_egenreco = new TH1F("h_egenreco","E reco sum over gen",100,0.,2.);//changed from 20 to 2
  TH1F* h_egensim = new TH1F("h_egensim","E sim sum over gen",100,0.,2.);//changed from 20 to 2

  TH2F* h_EpCone = new TH2F("h_EpCone","Ereco/gen versus cone size",10,0.,1.,100,0.,2.); // changed from 20 to 2 do to e weighting 
  TH2F* h_EpPhi = new TH2F("h_EpPhi","Ereco/gen versus phi",100,-4.,4.,100,0.,2.); // changed from 20 to 2 
  TH2F* h_etagenmax= new TH2F("h_etagenmax","eta gen vs max",100,1.,5.,100,1.,5.);
  TH2F* h_phigenmax= new TH2F("h_phigenmax","phi gen vs max",100,-4,4.,100,-4.,4.);
  TH1F* h_maxE = new TH1F("h_maxE","energy of highest energy hit",1000,0.,1000.);// changed from 5000 to 1000
  TH1F* h_ECone03 = new TH1F("h_ECone03","Sum energy cone 03",1000,0.,500.);// changed from 50000 to 500
  // histos from sarahs code //
  TH2F* h_banana = new TH2F("h_banana","banana plot",1000,0.,500.,1000,0.,500.);
  TH1F* h_fracBH = new TH1F("h_fracBH","fraction in BH",100,-01.,1.1);
  
  ///////////////////////////////////////////////////////
  //////////////////  start event loop
  //////////////////////////////////////////////////////

  //  const unsigned nEvts = ((pNevts > lRecTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lRecTree->GetEntries()) : pNevts) ;
  
  //std::cout << " -- Processing " << nEvts << " events out of " << 
  //   lRecTree->GetEntries()<< std::endl;


  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;

  std::cout << " -- Processing " << nEvts << " events out of " << lSimTree->GetEntries() << " " << lRecTree->GetEntries() << std::endl;

  //loop on events
  HGCSSEvent * event = 0;
  HGCSSEvent * eventRec = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lRecTree->SetBranchAddress("HGCSSEvent",&eventRec);
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

  //KH - myDigitiser
  HGCSSDigitisation myDigitiser;

  unsigned ievtRec = 0;
  unsigned nSkipped = 0;
  std::vector<double> absW; // from Sarah's code fill with weight
  bool firstEvent = true;

  // ---------- Event loop starts ----------

   for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  //for (unsigned ievt(0);ievt<100;++ievt){ // just lookings at the first 100 events for now
    if (ievtRec>=lRecTree->GetEntries()) continue;


    if (debug) std::cout << std::endl<<std::endl<<"... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;


    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievtRec);
    if (nPuVtx>0 && eventRec->eventNumber()==0 && event->eventNumber()!=0) {
      std::cout << " skip !" << ievt << " " << ievtRec << std::endl;
      nSkipped++;
      continue;
    }
    // start getting filling our weighting (from Sarah)

    if(firstEvent) {
      firstEvent=false;
      std::cout<<" size of ssvec of weights is "<<(*ssvec).size()<<std::endl;
      double absweight=0;      
      for (unsigned iL(0); iL<(*ssvec).size();++iL){
	if(iL<((*ssvec).size()-1)) {
	  unsigned next=iL+1;
	  absweight=(((*ssvec)[iL].voldEdx())+((*ssvec)[next].voldEdx()))/2. ;
	} else{
	  absweight+=(*ssvec)[iL].voldEdx()  ;
	}
	absW.push_back(absweight);
	absweight=0;
      }
      std::cout << " -- AbsWeight size: " << absW.size() << std::endl;
      std::cout<<" values are ";
      for (unsigned iL(0); iL<(*ssvec).size();++iL){
	//KH std::cout<<" "<<absW[iL];
	printf("layer# voldEdx, absW, volLambdatrans(), volX0trans(): %4d %8.3f %8.3f %8.3f %8.3f\n",
	       iL,(*ssvec)[iL].voldEdx(),absW[iL],(*ssvec)[iL].volLambdatrans(),(*ssvec)[iL].volX0trans());
      }
      std::cout<<std::endl;
    }
    // end of weightings

    double ptgen=-1.;
    double Egen=-1.;
    double ptgenpx=-1.;
    double ptgenpy=-1.;
    double ptgenpz=-1.;
    double etagen=99999.;
    double phigen=99999.;
    int pidgen=-1;
    double massgen= -1;
    if((*genvec).size()>0) {
      pidgen=(*genvec)[0].pdgid();
      massgen=(*genvec)[0].mass();
      ptgenpx=(*genvec)[0].px()/1000.;
      ptgenpy=(*genvec)[0].py()/1000.;
      ptgenpz=(*genvec)[0].pz()/1000.;
      ptgen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy);
      //double theta=atan(ptgen/ptgenpz);
      //etagen=-log(tan(theta/2));
      etagen=(*genvec)[0].eta();
      Egen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy+ptgenpz*ptgenpz);
      //phigen=atan2(ptgenpy,ptgenpx);
      phigen=(*genvec)[0].phi();
      //if(phigen<0) phigen=2.*TMath::Pi()+phigen;
    }
    h_getaphi->Fill(etagen,phigen);
    if(debug) {
      std::cout<<" gen vec size is "<<(*genvec).size()<<std::endl;
      std::cout<<" first gen   pt  "<<ptgen<<" egen  "<<Egen<<" pidgen  "<<pidgen<<" etagen  "<<etagen<<" phi gen "<<phigen<<  "mass gen "<< massgen<<std::endl;
      for (unsigned iP(0); iP<(*genvec).size(); ++iP){
        std::cout<<" gen particle "<<iP<<" is (pdgid) "<<(*genvec)[iP].pdgid()<<std::endl;
      }
    }


    bool isScint = false;
    if (debug) std::cout << " - Event contains " << (*rechitvec).size() << " rechits." << std::endl;

    // make some simple plots about all the rechits
    unsigned iMax=-1;
    double MaxE=-1.;

    //find rmin and rmax of rechits 
    double rmin = 99999;
    double rmax = -1;

    // ---------- Rechit loop starts ----------

    std::map<std::pair<int,int>,float> mymap_rechit;
    mymap_rechit.clear();

    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      double leta = lHit.eta();
      double lphi = lHit.phi();
      unsigned layer = lHit.layer();
      if (lHit.energy()>MaxE) {MaxE=lHit.energy(); iMax=iH;}
      if (debug>20) std::cout << " -- hit " << iH << " eta " << leta << std::endl; 

      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
      isScint = subdet.isScint;
      TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();

      unsigned cellid = map->FindBin(lHit.get_x(),lHit.get_y());
      geomConv.fill(lHit.layer(),lHit.energy(),0,cellid,lHit.get_z());
      
      if (debug>2 && lHit.energy()>1000.) std::cout << "reco energy"<< lHit.energy() << std::endl;
      //std::cout << "x "<< lHit.get_x() << "\t y "<<lHit.get_y() << "\t z" << lHit.get_z()<< std::endl; // added by Bryan, prints out xyz of each reco hit
      //std::cout<<"reco energy " << lHit.energy()<< " reco eta "<<lHit.eta()<< " reco phi " << lHit.phi() << " reco layer " << lHit.layer()<<" reco noise ratio "<< lHit.noiseFraction()<<std::endl;

      double lenergy=lHit.energy()*absW[layer]/1000.; // weight added (from Sarah's code)
      double r_hit = sqrt(lHit.get_x()*lHit.get_x()+lHit.get_y()*lHit.get_y());

      // printf added by bryan
      /*
      printf("|| reco hit# e, corre, eta, phi, lyr, noiseF: %5d %6.1f %6.1f %6.1f %6.1f %3d %8.3f\n",
	     iH,lHit.energy(),lenergy,lHit.eta(),lHit.phi(),lHit.layer(),lHit.noiseFraction());
      printf("|| reco hit# %d  \t",iH);
      printf("| reco energy = %f \t",lHit.energy());
      printf("| reco weighted E = %f \t", lenergy);
      printf("| reco eta = %f \t",lHit.eta());
      printf("| reco phi = %f \t",lHit.phi());
      printf("| reco layer = %d \t",lHit.layer());
      printf("| reco noise ratio = %f\t ||\n ",lHit.noiseFraction());
      */

      //if (z_layer[lHit.layer()]<1.)
      //	z_layer[lHit.layer()] = lHit.get_z();

      h_energy->Fill(lHit.energy());
      h_z->Fill(lHit.get_z());
      h_z1->Fill(lHit.get_z());
      h_z2->Fill(lHit.get_z());
      h_l->Fill(lHit.layer()+0.5);
      int ixx=lHit.layer();
      if(ixx>52) ixx=ixx-17;
      h_zl->Fill(lHit.get_z(),ixx);
      h_l2->Fill(lHit.layer()+0.5);
      h_xy->Fill(lHit.get_x(),lHit.get_y());
      h_zr->Fill(lHit.get_z(),r_hit); // added by Bryan
      h_zx->Fill(lHit.get_z(),lHit.get_x()); //added by Bryan
      h_yz->Fill(lHit.get_z(),lHit.get_y()); //added by Bryan
      h_zx1->Fill(lHit.get_z(),lHit.get_x()); //added by Bryan
      h_yz1->Fill(lHit.get_z(),lHit.get_y()); //added by Bryan
      h_zx2->Fill(lHit.get_z(),lHit.get_x()); //added by Bryan
      h_yz2->Fill(lHit.get_z(),lHit.get_y()); //added by Bryan
      h_zx10000->Fill(lHit.get_z(),lHit.get_x());//added by Bryan
      h_zx1000->Fill(lHit.get_z(),lHit.get_x());//added by Bryan
      h_xyz->Fill(lHit.get_x(),lHit.get_y(),lHit.get_z());// added by Bryan
      h_etaphi->Fill(lHit.eta(),lHit.phi());

      if(!isScint && ixx == 37){
	if(r_hit > rmax){
	  rmax = r_hit;
	};
      }

      if(!isScint) {
	map_1->Fill(lHit.get_x(),lHit.get_y());
	if (develop && ixx>=36) map_1_layer[ixx]->Fill(lHit.get_x(),lHit.get_y());
      }
      //if (ixx<=20) map_1_0->Fill(lHit.get_x(),lHit.get_y());
      //else         map_1_1->Fill(lHit.get_x(),lHit.get_y());

      if(isScint){ 
	if (ixx>=36&&ixx<=39)      map_TH2F_2[ixx-36]->Fill(lHit.phi(),r_hit);
	else if (ixx>=40&&ixx<=51) map_TH2F_3[ixx-40]->Fill(lHit.phi(),r_hit);
      }

      if(isScint)
	{
	  h_sxy->Fill(lHit.get_x(),lHit.get_y());
	  h_szx->Fill(lHit.get_z(),lHit.get_x());
	}
      else
	{
	  h_nsxy->Fill(lHit.get_x(),lHit.get_y());
	  h_nszx->Fill(lHit.get_z(),lHit.get_x());
	}
      
      //initialise calibration class
      const unsigned nSiLayers = 2;  // this is what I see in generation, but why?
      if (debug>2) std::cout << "KHKH: inFilePath,bypassR,nSiLayers: " << inFilePath<<" "<<bypassR<<" "<<nSiLayers << std::endl;
      HGCSSCalibration mycalib(inFilePath,bypassR,nSiLayers);
      mycalib.setVertex(eventRec->vtx_x(),eventRec->vtx_y(),eventRec->vtx_z());

      //KH noise tests
      /*
      if(isScint){
	printf("Scinti (E,r,lay,adc): %8.2f %8.2f %8d %8.2f\n",lHit.energy(),r_hit,layer,lenergy/mycalib.MeVToMip(layer,r_hit));
      } else {
	printf("Si     (E,r,lay,adc): %8.2f %8.2f %8d %8.2f\n",lHit.energy(),r_hit,layer,lenergy/mycalib.MeVToMip(layer,r_hit));
      }
      */

      int ilayer = ixx;

      // added in zx comps. Bryan//////

      if(ilayer==36) {
	if(isScint)
	  {
	    h_sxy36->Fill(lHit.get_x(),lHit.get_y());
	    h_szx36->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else 
	  { 
	    h_nsxy36->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx36->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==37) {
	if(isScint)
	  {
	    h_sxy37->Fill(lHit.get_x(),lHit.get_y());
	    h_szx37->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy37->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx37->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==38) {
	if(isScint)
	  {
	    h_sxy38->Fill(lHit.get_x(),lHit.get_y());
	    h_szx38->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy38->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx38->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==39) {
	if(isScint)
	  {
	    h_sxy39->Fill(lHit.get_x(),lHit.get_y());
	    h_szx39->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy39->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx39->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==40) {
	if(isScint)
	  {
	    h_sxy40->Fill(lHit.get_x(),lHit.get_y());
	    h_szx40->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy40->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx40->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==41) {
	if(isScint)
	  {
	    h_sxy41->Fill(lHit.get_x(),lHit.get_y());
	    h_szx41->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy41->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx41->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==42) {
	if(isScint)
	  {
	    h_sxy42->Fill(lHit.get_x(),lHit.get_y());
	    h_szx42->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy42->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx42->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==43) {
	if(isScint)
	  {
	    h_sxy43->Fill(lHit.get_x(),lHit.get_y());
	    h_szx43->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy43->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx43->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==44) {
	if(isScint)
	  {
	    h_sxy44->Fill(lHit.get_x(),lHit.get_y());
	    h_szx44->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy44->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx44->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==45) {
	if(isScint)
	  {
	    h_sxy45->Fill(lHit.get_x(),lHit.get_y());
	    h_szx45->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy45->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx45->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==46) {
	if(isScint)
	  {
	    h_sxy46->Fill(lHit.get_x(),lHit.get_y());
	    h_szx46->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy46->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx46->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==47) {
	if(isScint)
	  {
	    h_sxy47->Fill(lHit.get_x(),lHit.get_y());
	    h_szx47->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy47->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx47->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==48) {
	if(isScint)
	  {
	    h_sxy48->Fill(lHit.get_x(),lHit.get_y());
	    h_szx48->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy48->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx48->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==49) {
	if(isScint)
	  {
	    h_sxy49->Fill(lHit.get_x(),lHit.get_y());
	    h_szx49->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy49->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx49->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==50) {
	if(isScint)
	  {
	    h_sxy50->Fill(lHit.get_x(),lHit.get_y());
	    h_szx50->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy50->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx50->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };
      if(ilayer==51) {
	if(isScint)
	  {
	    h_sxy51->Fill(lHit.get_x(),lHit.get_y());
	    h_szx51->Fill(lHit.get_z(),lHit.get_x());
	    h_szxl->Fill(lHit.get_z(),lHit.get_x());
	    h_sxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
	else
	  {
	    h_nsxy51->Fill(lHit.get_x(),lHit.get_y());
	    h_nszx51->Fill(lHit.get_z(),lHit.get_x());
	    h_nszxl->Fill(lHit.get_z(),lHit.get_x());
	    h_nsxyl->Fill(lHit.get_x(),lHit.get_y());
	  }
      };


    }//loop on hits

    HGCSSRecoHit lHit = (*rechitvec)[iMax];
    double maxeta = lHit.eta();
    double maxphi=lHit.phi();
    double maxE=lHit.energy();
    h_maxE->Fill(maxE);
    h_etagenmax->Fill(maxeta,etagen);
    h_phigenmax->Fill(maxphi,phigen);
    if(debug>2) {
      std::cout<<" Max hit energy eta phi "<<lHit.energy()<<" "<<lHit.eta()<<" "<<lHit.phi()<<std::endl;
    }

    // make e/p plots for various cones around gen particle
    double rechitsumE01=0.;
    double rechitsumE02=0.;
    double rechitsumE03=0.;
    double rechitsumE04=0.;
    double rechitsumE05=0.;
    double rechitBHsumE01=0.;// if a rechit is part of the scint 
    double rechitBHsumE02=0.;
    double rechitBHsumE03=0.;
    double rechitBHsumE04=0.;
    double rechitBHsumE05=0.;
    double rechitsumEWoNoise01=0.;
    double rechitsumEWoNoise02=0.;
    double rechitsumEWoNoise03=0.;
    double rechitsumEWoNoise04=0.;
    double rechitsumEWoNoise05=0.;
    double etaW=0.;
    double phiW=0.;
    double norm=0.;
    //    double etaaxis=etagen;
    //    double phiaxis=phigen;
    double etaaxis=maxeta;
    double phiaxis=maxphi;
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      unsigned layer = lHit.layer();
      double leta = lHit.eta();
      double lphi = lHit.phi();
      //if(lphi<0) lphi=2.*TMath::Pi()+lphi;
      double lenergy=lHit.energy()*absW[layer]/1000.; // weight added (from Sarah's code)
      if (debug>20) std::cout << " -- hit " << iH << " et eta phi " << lenergy<<" "<<leta << " "<< lphi<<std::endl; 
	//clean up rechit collection

      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
      isScint = subdet.isScint;

      norm+=lenergy;
      etaW+=leta*lenergy;
      phiW+=lphi*lenergy;

      double dR=DeltaR(etaaxis,phiaxis,leta,lphi);
      //double dR=fabs(etagen-leta);

      TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();

      unsigned cellid = map->FindBin(lHit.get_x(),lHit.get_y());

      //
      //KH - make a map for selected rechits around gen pions and low noise fraction
      //KH - Warning - cellid not working well for scintilator layers (so not counted)
      //
      if (!isScint && lHit.noiseFraction()<0.5 && dR<0.3){
	  unsigned cellid = map->FindBin(lHit.get_x(),lHit.get_y());
	  std::map<std::pair<int,int>, float>::iterator it = mymap_rechit.find(std::make_pair(layer,cellid)); 
	  if (it != mymap_rechit.end()) (*it).second += lHit.energy();
	  else mymap_rechit.insert(std::make_pair(std::make_pair(layer,cellid),lHit.energy()));
      }

      if(debug>20) std::cout<<" dR "<<dR<<" "<<etagen<<" "<<phigen<<" "<<leta<<" "<<lphi<<std::endl;
      if(dR<0.1)
	{
	  rechitsumE01+=lenergy;
	  if (isScint) rechitBHsumE01+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise01+=lenergy;
	}
      if(dR<0.2)
	{
	  rechitsumE02+=lenergy;
	  if (isScint) rechitBHsumE02+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise02+=lenergy;
	}
      if(dR<0.3)
	{
	  rechitsumE03+=lenergy;
	  if (isScint) rechitBHsumE03+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise03+=lenergy;
	}
      if(dR<0.4)
	{
	  rechitsumE04+=lenergy;
	  if (isScint) rechitBHsumE04+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise04+=lenergy;
	}
      if(dR<0.5)
	{
	  rechitsumE05+=lenergy;
	  if (isScint) rechitBHsumE05+=lenergy;
	  if (lHit.noiseFraction()<0.5) rechitsumEWoNoise05+=lenergy;
	}

    }//loop on hits
    if(debug>1) {
      std::cout<<" reco gen are "<<rechitsumE01<<" "<<rechitsumE02<<" "<<rechitsumE03<<" "<<rechitsumE04<<" "<<rechitsumE05<<" "<<Egen<<std::endl;
    }
    if(debug>5) std::cout<<"weighted eta phi are "<<etaW/norm<<" "<<phiW/norm<<std::endl;

    if (MaxE>energy_max) energy_max=MaxE;

    h_Egenreco->Fill(Egen,rechitsumE03/Egen);// changed from 5 to 3
    h_egenreco->Fill(rechitsumE03/Egen); // changed from 5 to 3
    h_EpPhi->Fill(phigen,rechitsumE03/Egen); // changed from 2 to 3
    h_EpCone->Fill(0.1,rechitsumE01/Egen);
    h_EpCone->Fill(0.2,rechitsumE02/Egen);
    h_EpCone->Fill(0.3,rechitsumE03/Egen);
    h_EpCone->Fill(0.4,rechitsumE04/Egen);
    h_EpCone->Fill(0.5,rechitsumE05/Egen);
    
    h_banana->Fill(rechitsumE03-rechitBHsumE03,rechitBHsumE03); // added from sarah's code
    double frac =-0.05; // from sarah
    double notBH=rechitsumE03-rechitBHsumE03; // from sarah unused 
    if(rechitsumE03>0) frac=rechitBHsumE03/rechitsumE03; // from sarah
    h_fracBH->Fill(frac); // from sarah

    h_ECone03->Fill(rechitsumE03);
      //miptree->Fill();

    // ---------- Simhit loop starts ----------

    double simhitsumE01=0.;
    double simhitsumE02=0.;
    double simhitsumE03=0.;
    double simhitsumE04=0.;
    double simhitsumE05=0.;
    double simhitBHsumE01=0.;
    double simhitBHsumE02=0.;
    double simhitBHsumE03=0.;
    double simhitBHsumE04=0.;
    double simhitBHsumE05=0.;

    //initialise calibration class
    //KH - myDigitiser
    //HGCSSDigitisation myDigitiser; // Now outside the event loop
    //const unsigned interCalib = 3; // check against generation setting    
    //myDigitiser.setIntercalibrationFactor(interCalib);
    const unsigned nSiLayers = 2;  // this is what I see in generation, but why?
    if (debug>2) std::cout << "KHKH: inFilePath,bypassR,nSiLayers: " << inFilePath<<" "<<bypassR<<" "<<nSiLayers << std::endl;
    HGCSSCalibration mycalib(inFilePath,bypassR,nSiLayers);
    mycalib.setVertex(eventRec->vtx_x(),eventRec->vtx_y(),eventRec->vtx_z());

    int nhit=0, nhitscinti=0;
    std::map<std::pair<int,int>,float> mymap_simhit; mymap_simhit.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_x; mymap_simhit_x.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_y; mymap_simhit_y.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_z; mymap_simhit_z.clear();    
    //std::map<std::pair<int,int>,float> mymap_simhit_eta; mymap_simhit_eta.clear();
    //std::map<std::pair<int,int>,float> mymap_simhit_phi; mymap_simhit_phi.clear();
    std::map<std::pair<int,int>,float> mymap_simhit_MeVToMip; mymap_simhit_MeVToMip.clear();

    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];

      unsigned layer = lHit.layer();
      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
      DetectorEnum type = subdet.type;
      isScint = subdet.isScint;
    
      std::pair<double,double> xy = lHit.get_xy(subdet,geomConv,shape);
      double posx = xy.first;//lHit.get_x(cellSize);
      double posy = xy.second;//lHit.get_y(cellSize);
      double posz = lHit.get_z();
      double radius = sqrt(pow(posx,2)+pow(posy,2));
      double energy = lHit.energy()*mycalib.MeVToMip(layer,radius); // if (energy > 0) std::cout << "sim energy = "<<lHit.energy()<<", reco energy = "<<energy<<std::endl;
      double absweight = absW[layer];
      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);

      //KH - myDigitiser
      bool passTime = myDigitiser.passTimeCut(type,realtime);
      if (!passTime) continue;
      if (realtime<-10.) continue;
      if (lHit.time()==0.) continue; // remove bad timing hits
      // This is not a rigorous timing cut, but clean some simhits which won't make visible differences in results

      ROOT::Math::XYZPoint lpos = ROOT::Math::XYZPoint(posx,posy,posz);
      double eta = lpos.eta();
      double phi = lpos.phi();
      double dR=DeltaR(etaaxis,phiaxis,eta,phi);

      if (!isScint){ // si
	DetectorEnum adet = subdet.type;
	unsigned adc = 0;
	double digiE;
	adc = myDigitiser.adcConverter(energy,type);
	/* KH debug
	digiE = myDigitiser.adcToMIP(adc,type);	
	printf("Si: adc MIP %8d %8.2f\n",adc,digiE);
	*/
	if (adc<5) continue;                // ADC cut on silicon hits
      } else {
	/* KH debug 
	DetectorEnum adet = subdet.type;
	unsigned adc = 0;
	double digiE;
	adc = myDigitiser.adcConverter(energy,type);
	adc = myDigitiser.adcConverter(energy,type);
	digiE = myDigitiser.adcToMIP(adc,type);	
	printf("Scinti: adc MIP %8d %8.2f\n",adc,digiE);
	*/
	if (energy<0.5) continue;           // MIP cut (>=0.5MIP) on scintilator hits
      }

      // 
      if (layer<=51&&lHit.silayer()>=geomConv.getNumberOfSiLayers(type,radius)) continue;

      nhit++;
      if (layer>51) nhitscinti++;
      //-----
      /*
      if (layer>51)
      std::cout << "MeVToMip: " << layer << " "
		<< radius << " "
		<< lHit.cellid() << " "
		<< eta << ", " << phi << ", " << dR << " " << energy << " "
		<< mycalib.MeVToMip(layer,radius) << std::endl;
      */
      
      if (!isScint && dR<0.3){
	std::map<std::pair<int,int>, float>::iterator it = mymap_simhit.find(std::make_pair(layer,lHit.cellid())); 
	if (it != mymap_simhit.end()){
	  (*it).second += energy;
	}
	else {
	  mymap_simhit.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),energy));
	  mymap_simhit_MeVToMip.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),mycalib.MeVToMip(layer,radius)));
	  mymap_simhit_x.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),posx));
	  mymap_simhit_y.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),posy));
	  mymap_simhit_z.insert(std::make_pair(std::make_pair(layer,lHit.cellid()),posz));
	}
      }
      
      if(dR<0.1){ 
	simhitsumE01+=energy*absweight/1000.;
	if (isScint) simhitBHsumE01+=energy*absweight/1000.;
      }
      if(dR<0.2){
	simhitsumE02+=energy*absweight/1000.;
	if (isScint) simhitBHsumE02+=energy*absweight/1000.;
      }
      if(dR<0.3){
	simhitsumE03+=energy*absweight/1000.;
	if (isScint) simhitBHsumE03+=energy*absweight/1000.;
      }
      if(dR<0.4){
	simhitsumE04+=energy*absweight/1000.;
	if (isScint) simhitBHsumE04+=energy*absweight/1000.;
      }
      if(dR<0.5){
	simhitsumE05+=energy*absweight/1000.;      
	if (isScint) simhitBHsumE05+=energy*absweight/1000.;
      }
      //std::cout << absWeight << " " << absW[layer] << std::endl;

    }

    /////-----
    /*
    bool print_header=true;
    for (std::map<std::pair<int,int>,float>::iterator it=mymap_simhit.begin(); it!=mymap_simhit.end(); ++it){
      std::pair<int,int> key=it->first;

      std::map<std::pair<int,int>, float>::iterator iit = mymap_rechit.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator iiit = mymap_simhit_MeVToMip.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator ix = mymap_simhit_x.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator iy = mymap_simhit_y.find(std::make_pair(key.first,key.second)); 
      std::map<std::pair<int,int>, float>::iterator iz = mymap_simhit_z.find(std::make_pair(key.first,key.second)); 
      double rr = pow(pow(ix->second,2)+pow(iy->second,2),0.5);
      if (it->second>0.){
      if (print_header){
	printf("mymap_simhit: cellid, layer, sim E (mip) [rec E, rec/sim E, MeVToMip, rr]\n");
	print_header=false;
      }
      if (iit == mymap_rechit.end()) printf("mymap_simhit: %6d %10d %10.4f\n",key.first,key.second,it->second);
      else{ 
	printf("mymap_simhit: %6d %10d %10.4f %10.4f %10.4f %10.4f %8.3f\n",
	       key.first,key.second,
	       it->second,iit->second,iit->second/it->second,iiit->second,rr);      
	if (iit->second/it->second<0.8 && it->second>50.){
	  
	  for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	    HGCSSSimHit lHit = (*simhitvec)[iH];
	    if (lHit.layer()==key.first&&lHit.cellid()==key.second){

	      unsigned layer = lHit.layer();
	      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	      DetectorEnum type = subdet.type;
	      isScint = subdet.isScint;
	      std::pair<double,double> xy = lHit.get_xy(subdet,geomConv,shape);
	      double posx = xy.first;//lHit.get_x(cellSize);
	      double posy = xy.second;//lHit.get_y(cellSize);
	      double posz = lHit.get_z();
	      double radius = sqrt(pow(posx,2)+pow(posy,2));
	      double energy = lHit.energy()*mycalib.MeVToMip(layer,radius);
	      double absweight = absW[layer];
	      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
	      ROOT::Math::XYZPoint lpos = ROOT::Math::XYZPoint(posx,posy,posz);
	      double eta = lpos.eta();
	      double phi = lpos.phi();
	      double dR=DeltaR(etaaxis,phiaxis,eta,phi);
	      
	      std::cout << "Details: " << layer << " "
			<< radius << " "
			<< lHit.cellid() << " "
			<< eta << ", " << phi << ", " << dR << " " << energy << " "
			<< realtime << " "
			<< posx << " " << posy << " " << posz << " "
			<< mycalib.MeVToMip(layer,radius) << std::endl;

	    } // strange layer#, cell	    	   
	  }   // loop over simhits
	}     // reco/sim<0.8 and sim energy>50mips
      }       // found reco/sim match
      }

      //printf("mymap_simhit.insert: %6d %10d %10.4f\n",key.first,key.second,it->second);
    }
    std::cout << "mymap_simhit.size: " << mymap_simhit.size() << std::endl;
    /////-----
    print_header=true;
    for (std::map<std::pair<int,int>,float>::iterator it=mymap_rechit.begin(); it!=mymap_rechit.end(); ++it){
      std::pair<int,int> key=it->first;
      std::map<std::pair<int,int>, float>::iterator iit = mymap_simhit.find(std::make_pair(key.first,key.second)); 
      if (it->second>0.){
      if (print_header){
	printf("mymap_rechit: cellid, layer, rec E (mip) [sim E, sim/rec E, MeVToMip, rr]\n");
	print_header=false;
      }
      if (iit == mymap_simhit.end()) printf("mymap_rechit: %6d %10d %10.4f\n",key.first,key.second,it->second);
      else printf("mymap_rechit: %6d %10d %10.4f %10.4f %10.4f\n",key.first,key.second,it->second,iit->second,iit->second/it->second);
      }
    }
    std::cout << "mymap_rechit.size: " << mymap_rechit.size() << std::endl;
    //std::cout << "mymap_rechit: " << it->first << " => " << it->second << '\n';
    */

    h_egensim->Fill(simhitsumE03/Egen);
    
    if (debug>2){
    printf("simhitsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",simhitsumE01,simhitsumE02,simhitsumE03);
    printf("rechitsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",rechitsumE01,rechitsumE02,rechitsumE03);
    printf(" (w/o noise):       %8.3f, %8.3f, %8.3f\n",rechitsumE01,rechitsumE02,rechitsumE03);
    printf("rec/sim sumE01,2,3: %8.3f, %8.3f, %8.3f\n",rechitsumE01/simhitsumE01,rechitsumE02/simhitsumE02,
	   rechitsumE03/simhitsumE03);
    printf(" (w/o noise):       %8.3f, %8.3f, %8.3f\n",rechitsumEWoNoise01/simhitsumE01,
	   rechitsumEWoNoise02/simhitsumE02,
	   rechitsumEWoNoise03/simhitsumE03);
    printf("simhitBHsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",simhitBHsumE01,simhitBHsumE02,simhitBHsumE03);
    printf("rechitBHsumE01,2,3:   %8.3f, %8.3f, %8.3f\n",rechitBHsumE01,rechitBHsumE02,rechitBHsumE03);
    printf("rec/sim sumE01,2,3: %8.3f, %8.3f, %8.3f\n",rechitBHsumE01/simhitBHsumE01,
	   rechitBHsumE02/simhitBHsumE02,
	   rechitBHsumE03/simhitBHsumE03);
    std::cout << "nhit: "  << nhit << " / " << (*simhitvec).size() << " " << nhitscinti << std::endl;
   }

    //=========

    geomConv.initialiseHistos();
    ievtRec++;
   }//loop on entries
  // ---------- Event loop ends ----------

  std::cout<< "max energy "<< energy_max << std::endl;

  std::cout<<" size of gen hits "<< (*genvec).size()<< " size of rechits " << (*rechitvec).size()<<std::endl;

  if(debug) std::cout<<"writing files"<<std::endl;

  // draw th2poly on canvas

  if (develop){

    map_1->GetZaxis()->SetRangeUser(1,200);
    map_1->Print();
    map_1->SetMinimum(0.5);
    map_1->Write();
    
    for (int ilayer=36;ilayer<=51;ilayer++){
      //sprintf(title,"map_1_%d",ilayer);
      map_1_layer[ilayer]->Print();
      map_1_layer[ilayer]->Write();
    }
  }
  
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;
  
}//main
