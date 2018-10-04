#include "fwdet_res.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hvectorcand.h"
#include "hvectorcandsim.h"
#include "hparticlecandsim.h"
#include "hparticlecand.h"
#include "hparticletool.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <sstream>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";
#define PI TMath::Pi()

using namespace std;


HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);
HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & beamVector, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
  if (!loop->setInput(""))
    {                                                    // reading file structure
      std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  TStopwatch timer;
  timer.Reset();
  timer.Start();

  //////////////////////////////////////////////////////////////////////////////
  //      Fast tree builder for creating of ntuples                            //
  //////////////////////////////////////////////////////////////////////////////

  loop->printCategories();    // print all categories found in input + status

  HCategory * fCatGeantKine = nullptr;
  fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");
  if (!fCatGeantKine)
    {
      cout << "No catGeantKine!" << endl;
      exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

  HCategory * fFwDetStrawCal = nullptr;
  fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");
  if (!fFwDetStrawCal)
    {
      cout << "No catFwDetStrawCal!" << endl;
      exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

  HCategory * fCatVectorCandSim = nullptr;
  fCatVectorCandSim = HCategoryManager::getCategory(catVectorCand, kTRUE, "catVectorCand");
  if (!fCatVectorCandSim)
    {
      cout << "No catVectorCand!" << endl;
      //exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

  HCategory * fCatParticleCandSim= nullptr;
  fCatParticleCandSim = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
  if(!fCatParticleCandSim)
    {
      cout<< "No catParticleCandSim!"<<endl;
    }

  Int_t entries = loop->getEntries();
  //     //setting numbers of events regarding the input number of events by the user
  if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;

  //     // specify output file
  TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
  output_file->cd();
  //
  cout << "NEW ROOT TREE " << endl;
  //
  //crete histograms
  TCanvas* cAngles=new TCanvas("cAngles","Angles simulated and reconstruced");
  TH1F* hThetaSim=new TH1F("hThetaSim","Theta simulated in vertex position",1000,0,0.4);
  TH1F* hThetaRec=new TH1F("hThetaRec","Theta reconstructed in detector",1000,0,0.4);
  TH1F* hPhiSim=new TH1F("hPhiSim","Phi simulated in vertex position",1000,-PI*2,2*PI);
  TH1F* hPhiRec=new TH1F("hPhiRec","Phi reconstructed in detector",1000,-PI*2,2*PI);

  
  TCanvas* cRes=new TCanvas("cRes", "resolution");
  TH1F* hThetaRes=new TH1F("hThetaRes","Theta angle resolution",1000,-0.04,0.04);
  TH1F* hPhiRes=new TH1F("hPhiRes","Phi angle resolution",1000,-0.04,0.04);
  TH2F* h2PhiRecVsSim=new TH2F("h2PhiRecVsSim","Phi reconstructed vs. simulated",300,0,2*PI,300,0,2*PI);
  TH2F* h2ThetaRecVsSim=new TH2F("h2ThetaVsSim","Theta reconstructed vs. simulated",300,0,0.3,300,0,0.3);

  TCanvas* cResCuts=new TCanvas("cResCuts","Resolution for theta<0.1");
  TH1F* hThetaResCuts=new TH1F("hThetaResCuts","Theta angle resolution",1000,-0.04,0.04);
  TH1F* hPhiResCuts=new TH1F("hPhiResCuts","Phi angle resolution",1000,-0.04,0.04);
  TH2F* h2PhiRecVsSimCuts=new TH2F("h2PhiRecVsSimCuts","Phi reconstructed vs. simulated;sim ;rec",300,0,2*PI,300,0,2*PI);
  TH2F* h2ThetaRecVsSimCuts=new TH2F("h2ThetaVsSimCuts","Theta reconstructed vs. simulated; sim; rec",300,0,0.15,300,0,0.15);

  TCanvas* cResInTh=new TCanvas("cResInTh","Resolution in function of theta");
  const int nsec=10;
  TH1F* hThetaResInTh[nsec];
  TH1F* hPhiResInTh[nsec];
  TF1* fGausPhi[nsec];
  TF1* fGausTh[nsec];
  int resolution=200;
  double xmin=-0.4;
  double xmax=-xmin;
  for(int i =0;i<nsec;i++)
    {
      ostringstream histogramNamePhi;
      histogramNamePhi << "PhiRes_area_no_" << i;
      ostringstream histogramNameTh;
      histogramNameTh << "ThetaRes_area_no_" << i;
      ostringstream gausphi;
      gausphi << "phi_gaus_" << i;
      ostringstream gausth;
      gausth << "theta_gaus_" << i;
      hThetaResInTh[i]=new TH1F(histogramNameTh.str().c_str(),"",resolution*4,xmin/5,xmax/5);
      hPhiResInTh[i]=new TH1F(histogramNamePhi.str().c_str(),"",resolution,xmin,xmax);
      fGausPhi[i]=new TF1(gausphi.str().c_str(),"gaus",xmin,xmax);
      fGausTh[i]=new TF1(gausth.str().c_str(),"gaus",xmin/20,xmax/20);
    }
  
  TCanvas* cResXY=new TCanvas("cResXY", "Resolution for XY measurment");
  TH1F* hRX=new TH1F("hRX","Rsolution in X-direction for first plane",500,-5,5);
  TH1F* hRY=new TH1F("hRY","Rsolution in Y-direction for first plane",500,-5,5);
  TH2F* h2XRecVsSim=new TH2F("h2XRecVsSim","X reconstructed vs. simulated;sim;rec",200,-400,400,200,-400,400);
  TH2F* h2YRecVsSim=new TH2F("h2YRecVsSim","Y reconstructed vs. simulated;sim;rec",200,-400,400,200,-400,400);

  TCanvas* cXY=new TCanvas("cXY","X-Y hit distribution");
  TH1F* hXdistr=new TH1F("hXdistr","X hit distribution reconstructed",400,-400,400);
  TH1F* hYdistr=new TH1F("hYdistr","Y hit distribution reconstructed",400,-400,400);
  TH1F* hXsim=new TH1F("hXsim","X hit distribution simulated",400,-400,400);
  TH1F* hYsim=new TH1F("hYsim","Y hit distribution simulated",400,-400,400);

  TCanvas* cGeant=new TCanvas("cGeant","Points saved from GEANT");
  TH2F* h2PointX1=new TH2F("h2PointX1","First point registered in FwDet;z;x",200,3140,3240,2000,-500,500);
  TH2F* h2PointY1=new TH2F("h2PointY1","First point registered in FwDet;z;y",100,5278,5294,2000,-500,500);
  TH2F* h2PointX2=new TH2F("h2PointX2","Last point registered in FwDet;z;x",200,3140,3240,2000,-500,500);
  TH2F* h2PointY2=new TH2F("h2PointY2","Last point registered in FwDet;z;y",100,5278,5294,2000,-500,500);

  TCanvas* cResInZ=new TCanvas("cResInZ","Resoult in function Z");
  const int leyer=3;
  TH1F* hThetaResInZ[leyer];
  TF1* fThetaResInZ[leyer];
  double xxmin=-0.01;
  double xxmax=0.01;
  for(int i=0;i<leyer;i++)
    {
      ostringstream histogramNameTh;
      histogramNameTh << "hThetaResInZ_" << i;
      ostringstream gausth;
      gausth << "fThetaResInZ_" << i;
      hThetaResInZ[i]=new TH1F(histogramNameTh.str().c_str(),histogramNameTh.str().c_str(),500,xxmin,xxmax);
      fThetaResInZ[i]=new TF1(gausth.str().c_str(),"gaus",xxmin/10,xxmax/10);
    }

   
  TCanvas* cPlots=new TCanvas("cPlots","Different plots");
  TGraph* gDfOverDt=new TGraph(nsec);
  TGraph* gPredict=new TGraph(nsec);
  TGraph* gResInZ=new TGraph(leyer);
  TGraph* gResInTh=new TGraph(nsec);


  TCanvas *cResU=new TCanvas("cResU","Ucoordinate for every leyer");
  const int nl=16;
  TH1F* hUSingleLeyer[nl];
  for(int i=0;i<nl;i++)
    {
      ostringstream histogramNameTh;
      histogramNameTh << "hUcoordinate_leyer_" << i;
      //ostringstream gausth;
      //gausth << "hUcoordinate_leyer_" << i;
      hUSingleLeyer[i]=new TH1F(histogramNameTh.str().c_str(),histogramNameTh.str().c_str(),1200,-600,600);
      // fThetaResInZ[i]=new TF1(gausth.str().c_str(),"gaus",-400,400);
    }
  

  // TCanvas* cXY2D=new TCanvas("cXY2D","X and Y rec vr. sim");
  // TH2F* h2X2D=new TH2F("h2X2D","X rec vs. X sim",400,-400,400,400,-400,400);
  // TH2F* h2Y2D=new TH2F("h2Y2D","Y rec vs. Y sim",400,-400,400,400,-400,400);
  //event loop *************************************************
  //*********************************************************
  double thmin=0.03;
  double thmax=0.1;
  double thstep=(thmax-thmin)/nsec;
  for (Int_t i = 0; i < entries; i++)                   
    {
      loop->nextEvent(i);         // get next event. categories will be cleared before
      if(i%5000==0)
	cout<<"event no. "<<i<<endl;
      HParticleCandSim* particlecand =nullptr;
      HVectorCandSim* fwdetstrawvec = nullptr;
      HGeantKine* hkine=nullptr;
      HFwDetStrawCalSim* strawcal=nullptr;
      int vcnt, gknt;
      double phisim,phisim0,phisim1, phirec,thetasim, thetasim1,thetasim2,thetasim0, thetarec, xsim, ysim, xrec, yrec, zsts1, rsim,rsim1,rsim2;
   	
      //FW Detector
      if (fCatVectorCandSim)
	{
	  vcnt = fCatVectorCandSim->getEntries();
	  gknt = fCatGeantKine->getEntries();
	  for (int j = 0; j < vcnt; ++j)
	    {
	      fwdetstrawvec = HCategoryManager::getObject(fwdetstrawvec, fCatVectorCandSim, j);
	      Int_t vectorcandTID=fwdetstrawvec->getTrack();
	      Int_t track_number=0;
	      Int_t creationID=-1;
	      Int_t nHits=0;
	      
	      for(int k=0;k<gknt;k++)//find fitting vector form kine
		{
		  hkine = HCategoryManager::getObject(hkine, fCatGeantKine,k);
		  if(vectorcandTID==hkine->getTrack())
		    {
		      creationID=hkine->getMechanism(); 
		      // if(creationID !=0)//only primary vertex
		      //	continue;

		      rsim1=TMath::Sqrt(fwdetstrawvec->getPx1()*fwdetstrawvec->getPx1()+fwdetstrawvec->getPy1()*fwdetstrawvec->getPy1());
		      rsim2=TMath::Sqrt(fwdetstrawvec->getPx2()*fwdetstrawvec->getPx2()+fwdetstrawvec->getPy2()*fwdetstrawvec->getPy2());
		      rsim=rsim1;

		      thetasim0=hkine->getThetaDeg()*PI/180;
		      thetasim1=TMath::ATan(rsim1/fwdetstrawvec->getPz1());
		      thetasim2=TMath::ATan(rsim2/fwdetstrawvec->getPz2());
		      thetasim=thetasim1;

		      phisim1=-TMath::ATan(fwdetstrawvec->getPx1()/fwdetstrawvec->getPy1());
		      phisim0=hkine->getPhiDeg()*PI/180;
		      phisim=phisim0;
		      
		      thetarec=fwdetstrawvec->getHadesTheta();
		      phirec=fwdetstrawvec->getHadesPhi();

		      zsts1=fwdetstrawvec->getZ();
		      xrec=fwdetstrawvec->getX();
		      yrec=fwdetstrawvec->getY();
		      xsim=(zsts1-fwdetstrawvec->getZ1())*fwdetstrawvec->getPx1()/fwdetstrawvec->getPz1()+fwdetstrawvec->getX1();
		      ysim=(zsts1-fwdetstrawvec->getZ1())*fwdetstrawvec->getPy1()/fwdetstrawvec->getPz1()+fwdetstrawvec->getY1();

		      if(phirec<0)//redefine phi angle to range <0,2 Pi)
			phirec=2*PI+phirec;
		      //	      if(fwdetstrawvec->getPx1()<0)//redefine phi angle to range <0, Pi/2)
		      // phisim=phisim-PI;
		      // if(phisim<0)
		      // phisim=phisim+2*PI;

		      
		      hThetaSim->Fill(thetasim);
		      hPhiSim->Fill(phisim);
		      hThetaRec->Fill(thetarec);
		      hPhiRec->Fill(phirec);
		     	      
		      hThetaRes->Fill(thetarec-thetasim);
		      hPhiRes->Fill(phirec-phisim);
		      h2ThetaRecVsSim->Fill(thetasim,thetarec);
		      h2PhiRecVsSim->Fill(phisim,phirec);

		      h2XRecVsSim->Fill(xsim,xrec);
		      h2YRecVsSim->Fill(ysim,yrec);
		      hRY->Fill(yrec-ysim);
		      hRX->Fill(xsim-xrec);

		      hXdistr->Fill(xrec);
		      hYdistr->Fill(yrec);
		      hXsim->Fill(xsim);
		      hYsim->Fill(ysim);

		      h2PointY1->Fill(fwdetstrawvec->getZ1(),fwdetstrawvec->getY1());
		      h2PointY2->Fill(fwdetstrawvec->getZ2(),fwdetstrawvec->getY2());
		      h2PointX1->Fill(fwdetstrawvec->getZ1(),fwdetstrawvec->getX1());
		      h2PointX2->Fill(fwdetstrawvec->getZ2(),fwdetstrawvec->getX2());

		      //resoltion for every next layer
		      hThetaResInZ[0]->Fill(thetarec-thetasim0);
		      hThetaResInZ[1]->Fill(thetarec-thetasim1);
		      hThetaResInZ[2]->Fill(thetarec-thetasim2);

		      if(thetarec<0.1)//only phisical tracks
			{
			  hThetaResCuts->Fill(thetarec-thetasim);
			  hPhiResCuts->Fill(phirec-phisim);

			  h2ThetaRecVsSimCuts->Fill(thetasim,thetarec);
			  h2PhiRecVsSimCuts->Fill(phisim,phirec);
			}

		      //theta divided to 10 different crcles
		      for(int i=0;i<nsec;i++)
			{
			  if(thetasim>thmin+i*thstep && thetasim<=thmin+(i+1)*thstep)
			    {
			      hThetaResInTh[i]->Fill(thetarec-thetasim);
			      hPhiResInTh[i]->Fill(phirec-phisim);
			    }
			}

		      //residuals plot
		      nHits=fwdetstrawvec->getNofHits();
		      for(int i=0; i<nHits;i++)
			{
			  int layernumber;
			  int index=fwdetstrawvec->getHitIndex(i);
			  strawcal = HCategoryManager::getObject(strawcal,fFwDetStrawCal,index);
			  layernumber=8*strawcal->getStation()+2*strawcal->getLayer()+strawcal->getPlane();
			  hUSingleLeyer[layernumber]->Fill(strawcal->getU());			  
			}
		      
		      break;
		    }
		}
	     
	    }
	}
          
      //Kine
      /*
      if(fCatGeantKine)
	{
	  gknt=fCatGeantKine->getEntries();
	  for(int i=0;i<gknt;i++)
	    {
	      hkine = HCategoryManager::getObject(hkine, fCatGeantKine,i);
	    }
	}
      */
	
	
    } // end eventloop
  //***********************************************************************************

  //drawing histograms
  cAngles->Divide(2,2);
  cAngles->cd(1);
  hThetaRec->Draw();
  cAngles->cd(2);
  hThetaSim->Draw();
  cAngles->cd(3);
  hPhiRec->Draw();
  cAngles->cd(4);
  hPhiSim->Draw();

  cRes->Divide(2,2);
  cRes->cd(1);
  hThetaRes->Draw();
  cRes->cd(2);
  hPhiRes->Draw();
  cRes->cd(3);
  h2ThetaRecVsSim->Draw("COLZ");
  cRes->cd(4);
  h2PhiRecVsSim->Draw("COLZ");

  cResCuts->Divide(2,2);
  cResCuts->cd(1);
  hThetaResCuts->Draw();
  cResCuts->cd(2);
  hPhiResCuts->Draw();
  cResCuts->cd(3);
  h2ThetaRecVsSimCuts->Draw("COLZ");
  cResCuts->cd(4);
  h2PhiRecVsSimCuts->Draw("COLZ");

  cResInTh->Divide(5,4);
  for(int i=0; i<nsec;i++)
    {
      double xtheta=thmin+(i+0.5)*thstep;
      double xz=1575;
      double xr=xz*TMath::Tan(xtheta);
      if(i+nsec+1>20)
	break;
      cResInTh->cd(i+1);
      hPhiResInTh[i]->Draw();
      hPhiResInTh[i]->Fit(fGausPhi[i]);

      cResInTh->cd(i+nsec+1);
      fGausTh[i]->SetParameter(0,hThetaResInTh[i]->GetMaximum());
      fGausTh[i]->SetParameter(1,0);
      fGausTh[i]->SetParameter(2,0.0005);
      hThetaResInTh[i]->Draw();
      hThetaResInTh[i]->Fit(fGausTh[i],"R");

      gDfOverDt->SetPoint(i,xtheta,fGausPhi[i]->GetParameter(2)/fGausTh[i]->GetParameter(2));
      gPredict->SetPoint(i,xtheta,(xr*xr+xz*xz)/(xr*xz));
      gResInTh->SetPoint(i,xtheta,fGausTh[i]->GetParameter(2));
    }

  cResXY->Divide(2,2);
  cResXY->cd(1);
  hRY->Draw();
  cResXY->cd(2);
  hRX->Draw();
  cResXY->cd(3);
  h2YRecVsSim->Draw("COLZ");
  cResXY->cd(4);
  h2XRecVsSim->Draw("COLZ");

  cXY->Divide(2,2);
  cXY->cd(1);
  hXsim->Draw();
  cXY->cd(2);
  hYsim->Draw();
  cXY->cd(3);
  hXdistr->Draw();
  cXY->cd(4);
  hYdistr->Draw();

  cGeant->Divide(2,2);
  cGeant->cd(1);
  h2PointX1->Draw("COLZ");
  cGeant->cd(2);
  h2PointX2->Draw("COLZ");
  cGeant->cd(3);
  h2PointY1->Draw("COLZ");
  cGeant->cd(4);
  h2PointY2->Draw("COLZ");

  cResInZ->Divide(3);
  for(int i=0; i<leyer; i++)
    {
      cResInZ->cd(i+1);
      hThetaResInTh[i]->Draw();
      hThetaResInTh[i]->Fit(fThetaResInZ[i],"R");
      hThetaResInTh[i]->SetAxisRange(-0.02,0.02,"X");
      gResInZ->SetPoint(i,i,fThetaResInZ[i]->GetParameter(2));
    }

  cResU->Divide(4,4);
  for(int k=0; k<16; k++)
    {
      cResU->cd(k+1);
      hUSingleLeyer[k]->Draw();
    }


  cPlots->Divide(2,2);
  cPlots->cd(1);
  // create a 2-d histogram to define the range
  TH2F *hr = new TH2F("hr","Several graphs in the same pad",200,0,0.15,100,5,40);
  hr->SetXTitle("Theta (rad)");
  hr->SetYTitle("d phi / d theta");
  hr->Draw();  
  gDfOverDt->Draw(/*"AC*"*/"LP");
  gDfOverDt->SetMarkerColor(kRed);
  gDfOverDt->SetMarkerStyle(20);
  gPredict->Draw("LP");
  gPredict->SetMarkerStyle(21);
  gPredict->SetMarkerColor(kBlue);

  cPlots->cd(2);
  // create a 2-d histogram to define the range
  TH2F *hz = new TH2F("hz","Resolution in function of leyer;no ler;sigma[rad]",40,-1,7,50,0,0.005);
  hz->Draw();  
  gResInZ->Draw("LP");
  gResInZ->SetMarkerColor(kGreen);
  gResInZ->SetMarkerStyle(21);

  cPlots->cd(3);
  TH2F *hTh = new TH2F("hTh","Resolution in function of theta;theta[rad];sigma theta[rad]",200,0,0.15,100,0,0.002);
  hTh->Draw();  
  gResInTh->Draw("LP");
  gResInTh->SetMarkerColor(kGreen);
  gResInTh->SetMarkerStyle(21);


  
  //save histograms
  cAngles->Write();
  cRes->Write();
  cResCuts->Write();
  cResInTh->Write();
  cPlots->Write();
  cResXY->Write();
  cXY->Write();
  cGeant->Write();
  cResInZ->Write();
  cResU->Write();
  
  output_file->Close();
  cout << "writing root tree done" << endl;

  timer.Stop();
  timer.Print();

  return 0;
}
