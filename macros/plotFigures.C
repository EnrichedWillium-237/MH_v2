# include "TCanvas.h"
# include "TDirectory.h"
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
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const int nbinsETA = 16;
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0, -2.4,  0.0, -2.4,  0.4};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  0.0,  2.4, -0.4,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24",
                   "-24_0", "0_24", "-24_-4", "4_24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 16;
TString ANAL[] = {"N1MCm22SUB2", "N1MCp22SUB2", "N1MCm22SUB3", "N1MCp22SUB3",
                  "N1SUB2",      "N1ASUB2",     "N1BSUB2",     "N112SUB2",    "N112ASUB2",   "N112BSUB2",
                  "N1SUB3",      "N1ASUB3",     "N1BSUB3",     "N112SUB3",    "N112ASUB3",   "N112BSUB3"
};

TH1D * N1MCm22SUB2_pT[nbinsETA][ncentbins];
TH1D * N1MCp22SUB2_pT[nbinsETA][ncentbins];
TH1D * N1MCm22SUB3_pT[nbinsETA][ncentbins];
TH1D * N1MCp22SUB3_pT[nbinsETA][ncentbins];
TH1D * N1SUB2_pT[nbinsETA][ncentbins];
TH1D * N1ASUB2_pT[nbinsETA][ncentbins];
TH1D * N1BSUB2_pT[nbinsETA][ncentbins];
TH1D * N112SUB2_pT[nbinsETA][ncentbins];
TH1D * N112ASUB2_pT[nbinsETA][ncentbins];
TH1D * N112BSUB2_pT[nbinsETA][ncentbins];
TH1D * N1SUB3_pT[nbinsETA][ncentbins];
TH1D * N1ASUB3_pT[nbinsETA][ncentbins];
TH1D * N1BSUB3_pT[nbinsETA][ncentbins];
TH1D * N112SUB3_pT[nbinsETA][ncentbins];
TH1D * N112ASUB3_pT[nbinsETA][ncentbins];
TH1D * N112BSUB3_pT[nbinsETA][ncentbins];

TH1D * N1MCm22SUB2_eta[ncentbins];
TH1D * N1MCp22SUB2_eta[ncentbins];
TH1D * N1MCm22SUB3_eta[ncentbins];
TH1D * N1MCp22SUB3_eta[ncentbins];
TH1D * N1SUB2_eta[ncentbins];
TH1D * N1ASUB2_eta[ncentbins];
TH1D * N1BSUB2_eta[ncentbins];
TH1D * N112SUB2_eta[ncentbins];
TH1D * N112ASUB2_eta[ncentbins];
TH1D * N112BSUB2_eta[ncentbins];
TH1D * N1SUB3_eta[ncentbins];
TH1D * N1ASUB3_eta[ncentbins];
TH1D * N1BSUB3_eta[ncentbins];
TH1D * N112SUB3_eta[ncentbins];
TH1D * N112ASUB3_eta[ncentbins];
TH1D * N112BSUB3_eta[ncentbins];

TFile * tfinPt;
TFile * tfinEta;

void plotFigures()
{
    tfinPt = new TFile("results/MH_combined_Pt.root","read");

    for (int ebin = 0; ebin<netabins; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString mtag = Form("MH_nominal/eta_%s/%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);

            N1MCm22SUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1MCm22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1MCp22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCm22SUB3_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1MCm22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB3_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1MCp22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N112SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N112A2SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB2_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N112BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB3_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB3_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB3_pT[ebin][cbin] = (TH1D *) tfinPt->Get(Form("%s/N1BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
        }
    }
    

}
