# include "TCanvas.h"
# include "TFile.h"
# include "TGraphErrors.h"
# include "TH1.h"
# include "TH2.h"
# include "TLegend.h"
# include "TMath.h"
# include "TPaveText.h"
# include "Riostream.h"
# include <iostream>

using namespace std;

# include "src/HiEvtPlaneList.h"
# include "src/style.h"

using namespace hi;

static const int ncentbins = 13;
static const int cminCENT[] = {0,  5, 10, 15, 20, 25, 30, 35, 40, 50, 60,  0, 20,  60};
static const int cmaxCENT[] = {5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 20, 60, 100};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCm10SUB2", "N1MCm06SUB2", "N1MCm02SUB2",
                 "N1MCp02SUB2", "N1MCp06SUB2", "N1MCp10SUB2", "N1MCp14SUB2", "N1MCp18SUB2", "N1MCp22SUB2",
                 "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCm10SUB3", "N1MCm06SUB3", "N1MCm02SUB3",
                 "N1MCp02SUB3", "N1MCp06SUB3", "N1MCp10SUB3", "N1MCp14SUB3", "N1MCp18SUB3", "N1MCp22SUB3"
};

TH1D * gN1MCm22SUB2_pt[ncentbins];
TH1D * gN1MCm18SUB2_pt[ncentbins];
TH1D * gN1MCm14SUB2_pt[ncentbins];
TH1D * gN1MCm10SUB2_pt[ncentbins];
TH1D * gN1MCm06SUB2_pt[ncentbins];
TH1D * gN1MCm02SUB2_pt[ncentbins];

TH1D * gN1MCp22SUB2_pt[ncentbins];
TH1D * gN1MCp18SUB2_pt[ncentbins];
TH1D * gN1MCp14SUB2_pt[ncentbins];
TH1D * gN1MCp10SUB2_pt[ncentbins];
TH1D * gN1MCp06SUB2_pt[ncentbins];
TH1D * gN1MCp02SUB2_pt[ncentbins];

TH1D * gN1MC22SUB2_pt[ncentbins];
TH1D * gN1MC18SUB2_pt[ncentbins];
TH1D * gN1MC14SUB2_pt[ncentbins];
TH1D * gN1MC10SUB2_pt[ncentbins];
TH1D * gN1MC06SUB2_pt[ncentbins];
TH1D * gN1MC02SUB2_pt[ncentbins];

TH1D * gN1MCm22SUB3_pt[ncentbins];
TH1D * gN1MCm18SUB3_pt[ncentbins];
TH1D * gN1MCm14SUB3_pt[ncentbins];
TH1D * gN1MCm10SUB3_pt[ncentbins];
TH1D * gN1MCm06SUB3_pt[ncentbins];
TH1D * gN1MCm02SUB3_pt[ncentbins];

TH1D * gN1MCp22SUB3_pt[ncentbins];
TH1D * gN1MCp18SUB3_pt[ncentbins];
TH1D * gN1MCp14SUB3_pt[ncentbins];
TH1D * gN1MCp10SUB3_pt[ncentbins];
TH1D * gN1MCp06SUB3_pt[ncentbins];
TH1D * gN1MCp02SUB3_pt[ncentbins];

TH1D * gN1MC22SUB3_pt[ncentbins];
TH1D * gN1MC18SUB3_pt[ncentbins];
TH1D * gN1MC14SUB3_pt[ncentbins];
TH1D * gN1MC10SUB3_pt[ncentbins];
TH1D * gN1MC06SUB3_pt[ncentbins];
TH1D * gN1MC02SUB3_pt[ncentbins];

void plotPtDistEven()
{
    
    TH1::SetDefaultSumw2();

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        gN1MCm22SUB2_eta[cbin] = new TH1D(Form("N1MCm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm18SUB2_eta[cbin] = new TH1D(Form("N1MCm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm14SUB2_eta[cbin] = new TH1D(Form("N1MCm14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm10SUB2_eta[cbin] = new TH1D(Form("N1MCm10SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm06SUB2_eta[cbin] = new TH1D(Form("N1MCm06SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm02SUB2_eta[cbin] = new TH1D(Form("N1MCm02SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1MCp22SUB2_eta[cbin] = new TH1D(Form("N1MCp22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp18SUB2_eta[cbin] = new TH1D(Form("N1MCp18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp14SUB2_eta[cbin] = new TH1D(Form("N1MCp14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp10SUB2_eta[cbin] = new TH1D(Form("N1MCp10SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp06SUB2_eta[cbin] = new TH1D(Form("N1MCp06SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp02SUB2_eta[cbin] = new TH1D(Form("N1MCp02SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1MC22SUB2_eta[cbin] = new TH1D(Form("N1MC22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC18SUB2_eta[cbin] = new TH1D(Form("N1MC18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC14SUB2_eta[cbin] = new TH1D(Form("N1MC14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC10SUB2_eta[cbin] = new TH1D(Form("N1MC10SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC06SUB2_eta[cbin] = new TH1D(Form("N1MC06SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC02SUB2_eta[cbin] = new TH1D(Form("N1MC02SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1MCm22SUB3_eta[cbin] = new TH1D(Form("N1MCm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm18SUB3_eta[cbin] = new TH1D(Form("N1MCm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm14SUB3_eta[cbin] = new TH1D(Form("N1MCm14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm10SUB3_eta[cbin] = new TH1D(Form("N1MCm10SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm06SUB3_eta[cbin] = new TH1D(Form("N1MCm06SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCm02SUB3_eta[cbin] = new TH1D(Form("N1MCm02SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1MCp22SUB3_eta[cbin] = new TH1D(Form("N1MCp22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp18SUB3_eta[cbin] = new TH1D(Form("N1MCp18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp14SUB3_eta[cbin] = new TH1D(Form("N1MCp14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp10SUB3_eta[cbin] = new TH1D(Form("N1MCp10SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp06SUB3_eta[cbin] = new TH1D(Form("N1MCp06SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MCp02SUB3_eta[cbin] = new TH1D(Form("N1MCp02SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1MC22SUB3_eta[cbin] = new TH1D(Form("N1MC22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC18SUB3_eta[cbin] = new TH1D(Form("N1MC18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC14SUB3_eta[cbin] = new TH1D(Form("N1MC14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC10SUB3_eta[cbin] = new TH1D(Form("N1MC10SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC06SUB3_eta[cbin] = new TH1D(Form("N1MC06SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1MC02SUB3_eta[cbin] = new TH1D(Form("N1MC02SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
    }

}
