#include "fwdet_tests.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hvectorcand.h"

#include <TCanvas.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

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

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    //////////////////////////////////////////////////////////////////////////////
    //      Fast tree builder for creating of ntuples                            //
    //////////////////////////////////////////////////////////////////////////////

    loop->printCategories();    // print all categories found in input + status
    //     loop->printChain();            // print all files in the chain
    //     loop->Print();

    //     HEventHeader * header = loop->getEventHeader();
    HCategory * fCatGeantFwDet = nullptr;
    fCatGeantFwDet = HCategoryManager::getCategory(catFwDetGeantRaw, kTRUE, "catGeantFwDet");

    if (!fCatGeantFwDet)
    {
        cout << "No catGeantFwDet!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * fFwDetStrawCal = nullptr;
    fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");

    if (!fFwDetStrawCal)
    {
        cout << "No catFwDetStrawCal!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * fVectorCand = nullptr;
    fVectorCand = HCategoryManager::getCategory(catVectorCand, kTRUE, "catVectorCand");

    if (!fVectorCand)
    {
        cout << "No catVectorCand!" << endl;
//         exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    //
    Int_t entries = loop->getEntries();
    //     //setting numbers of events regarding the input number of events by the user
    if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;

    //     // specify output file
    TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
    output_file->cd();
    //
    cout << "NEW ROOT TREE " << endl;
    //

    const size_t mods = 2;
    const size_t layers = 4;

    TH2I * h_gmodXcellY_xy[mods][layers];
    TCanvas * c_gmodX[2];

    TH2I * h_dmodXcellY_straw_x[mods][layers];
    TCanvas * c_dmodX[2];

    TH2I * h_t_xy, * h_t_txty;
    TCanvas * c_t;

    char buff1[200];
    char buff2[200];
    for (uint i = 0; i < mods; ++i)
    {
        for (uint j = 0; j < layers; ++j)
        {
            sprintf(buff1, "h_gmod%dcell%d_xy", i, j);
            sprintf(buff2, "x-y correlation, mod=%d, cell=%d;x [mm];y [mm]", i, j);
            h_gmodXcellY_xy[i][j] = new TH2I(buff1, buff2, 60, -6, 6, 130, -650, 650);

            h_gmodXcellY_xy[i][j]->GetXaxis()->SetLabelSize(0.06);
            h_gmodXcellY_xy[i][j]->GetYaxis()->SetLabelSize(0.06);

            h_gmodXcellY_xy[i][j]->GetXaxis()->SetTitleSize(0.06);
            h_gmodXcellY_xy[i][j]->GetYaxis()->SetTitleSize(0.06);

            h_gmodXcellY_xy[i][j]->GetXaxis()->SetTitleOffset(1.05);
            h_gmodXcellY_xy[i][j]->GetYaxis()->SetTitleOffset(1.75);

            sprintf(buff1, "h_dmod%dcell%d_straw_x", i, j);
            sprintf(buff2, "x-straw correlation, mod=%d, cell=%d;x [mm];straw", i, j);
            h_dmodXcellY_straw_x[i][j] = new TH2I(buff1, buff2, 65, -650, 650, 60, 0, 120);

            h_dmodXcellY_straw_x[i][j]->GetXaxis()->SetLabelSize(0.06);
            h_dmodXcellY_straw_x[i][j]->GetYaxis()->SetLabelSize(0.06);

            h_dmodXcellY_straw_x[i][j]->GetXaxis()->SetTitleSize(0.06);
            h_dmodXcellY_straw_x[i][j]->GetYaxis()->SetTitleSize(0.06);

            h_dmodXcellY_straw_x[i][j]->GetXaxis()->SetTitleOffset(1.05);
            h_dmodXcellY_straw_x[i][j]->GetYaxis()->SetTitleOffset(1.75);
        }

        sprintf(buff1, "c_gmod%d", i);
        c_gmodX[i] = new TCanvas(buff1, buff1, 800, 800);
        c_gmodX[i]->DivideSquare(4);

        sprintf(buff1, "c_dmod%d", i);
        c_dmodX[i] = new TCanvas(buff1, buff1, 800, 800);
        c_dmodX[i]->DivideSquare(4);
    }

    h_t_xy = new TH2I("h_t_xy", "X-Y", 65, -650., 650., 65, -650., 650.);
    h_t_txty = new TH2I("h_t_txty", "Tx-Ty", 80, -0.2, 0.2, 80, -0.2, 0.2);

    c_t = new TCanvas("c_t", "c_t", 800, 400);
    c_t->DivideSquare(2);

    for (Int_t i = 0; i < entries; i++)                    // event loop
    {
        /*Int_t nbytes =*/  loop->nextEvent(i);         // get next event. categories will be cleared before
        //cout << fCatGeantFwDet->getEntries() << endl;

        std::vector<float> xcord[mods][layers];

        int geant_fwdet_cnt = fCatGeantFwDet->getEntries();

        HGeantFwDet * gfwdet = new HGeantFwDet;

        // kine
        for (int j = 0; j < geant_fwdet_cnt; ++j)
        {
            gfwdet = HCategoryManager::getObject(gfwdet, fCatGeantFwDet, j);

            Char_t mod, layer;
            Int_t geantCell;
            Int_t plane, cell;
            Float_t hx, hy, hz, px, py, pz, tof, len, E;

            gfwdet->getAddress(mod, layer, geantCell);
            gfwdet->getHit(hx, hy, hz, px, py, pz, tof, len, E);

            Int_t l_panel = (Int_t) geantCell/100 - 1;
            Int_t l_block = (Int_t) (geantCell%100)/10 - 1;
            Int_t l_straw = (Int_t) geantCell%10 - 1;

            if (l_block % 2 == 0)
                plane = 0; // TODO add module and layer dep
            else
                plane = 1; // TODO as above

            Int_t n_blocks = ( mod == 0 ? 5 : 7);
            Int_t n_straws = 8;

            // for a single panel, blocks/2 blocks for a plane =>  n_blocks >> 1
            // 
            //     | no of straws in a plane of panel | * panel number +
            //     | number of straw in current panel |
            cell = ((n_blocks >> 1) * n_straws ) * l_panel +
                (l_block >> 1) * n_straws + l_straw;

            h_gmodXcellY_xy[(int)mod][(int)layer]->Fill(hx, hz);
            xcord[(int)mod][(int)layer].push_back(hx);
        }

        // digitizer
        int fwdet_calsim_cnt = fFwDetStrawCal->getEntries();

        HFwDetStrawCalSim * fwdetstrawcs = new HFwDetStrawCalSim;

        for (int j = 0; j < fwdet_calsim_cnt; ++j)
        {
            fwdetstrawcs = HCategoryManager::getObject(fwdetstrawcs, fFwDetStrawCal, j);

            Char_t mod, l, p;
            Int_t straw;
            Float_t time, adc, radi, X, Z;
            Int_t StrawN;

            fwdetstrawcs->getAddress(mod, l, p, straw);
            fwdetstrawcs->getHit(time, adc, X, Z, StrawN);

            for (uint k = 0; k < xcord[(int)mod][(int)l].size(); ++k)
            {
                h_dmodXcellY_straw_x[(int)mod][(int)l]->Fill(xcord[(int)mod][(int)l][k], StrawN);
            }
        }

        if (fVectorCand)
        {
            Int_t vcnt = fVectorCand->getEntries();

            HVectorCand * fwdetstrawvec = new HVectorCand;

            for (int j = 0; j < vcnt; ++j)
            {
                fwdetstrawvec = HCategoryManager::getObject(fwdetstrawvec, fVectorCand, j);
    //             printf("Nof=%d\n", fwdetstrawvec->getNofHits());
    //             printf("x=%f y=%f tx=%f ty=%f\n",
    //                    fwdetstrawvec->getX(),
    //                    fwdetstrawvec->getY(),
    //                    fwdetstrawvec->getTx(),
    //                    fwdetstrawvec->getTy());

                h_t_xy->Fill(fwdetstrawvec->getX(), fwdetstrawvec->getY());
                h_t_txty->Fill(fwdetstrawvec->getTx(), fwdetstrawvec->getTy());
            }
        }
    } // end eventloop

    // Write
    // geant
    for (uint i = 0; i < mods; ++i)
    {
        for (uint j = 0; j < layers; ++j)
        {
            c_gmodX[i]->cd(1+j);
            gPad->SetLeftMargin(0.2);
            gPad->SetBottomMargin(0.2);

            h_gmodXcellY_xy[i][j]->Draw("colz");
            h_gmodXcellY_xy[i][j]->Write();
        }
        c_gmodX[i]->Write();
    }

    // digitizer
    for (uint i = 0; i < mods; ++i)
    {
        for (uint j = 0; j < layers; ++j)
        {
            // digitizer
            c_dmodX[i]->cd(1+j);
            gPad->SetLeftMargin(0.2);
            gPad->SetBottomMargin(0.2);

            h_dmodXcellY_straw_x[i][j]->Draw("colz");
            h_dmodXcellY_straw_x[i][j]->Write();
        }
        c_dmodX[i]->Write();
    }

    c_t->cd(1);
    h_t_xy->Draw("colz");
    c_t->cd(2);
    h_t_txty->Draw("colz");

    h_t_xy->Write();
    h_t_txty->Write();
    c_t->Write();

    output_file->Close();
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();

    return 0;
}
