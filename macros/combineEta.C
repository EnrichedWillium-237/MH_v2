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
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24", "-24_0", "0_-24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCp22SUB2", "N1MCp18SUB2", "N1MCp14SUB2",
                  "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCp22SUB3", "N1MCp18SUB3", "N1MCp14SUB3",
                  "N1SUB2",      "N1ASUB2",     "N1BSUB2",     "N112SUB2",    "N112ASUB2",   "N112BSUB2",
                  "N1SUB3",      "N1ASUB3",     "N1BSUB3",     "N112SUB3",    "N112ASUB3",   "N112BSUB3"
};

TGraphErrors * gN1MCm22SUB2[ncentbins];
TGraphErrors * gN1MCm18SUB2[ncentbins];
TGraphErrors * gN1MCm14SUB2[ncentbins];
TGraphErrors * gN1MCp22SUB2[ncentbins];
TGraphErrors * gN1MCp18SUB2[ncentbins];
TGraphErrors * gN1MCp14SUB2[ncentbins];
TGraphErrors * gN1SUB2[ncentbins];
TGraphErrors * gN1ASUB2[ncentbins];
TGraphErrors * gN1BSUB2[ncentbins];
TGraphErrors * gN112SUB2[ncentbins];
TGraphErrors * gN112ASUB2[ncentbins];
TGraphErrors * gN112BSUB2[ncentbins];

TGraphErrors * gN1MCm22SUB3[ncentbins];
TGraphErrors * gN1MCm18SUB3[ncentbins];
TGraphErrors * gN1MCm14SUB3[ncentbins];
TGraphErrors * gN1MCp22SUB3[ncentbins];
TGraphErrors * gN1MCp18SUB3[ncentbins];
TGraphErrors * gN1MCp14SUB3[ncentbins];
TGraphErrors * gN1SUB3[ncentbins];
TGraphErrors * gN1ASUB3[ncentbins];
TGraphErrors * gN1BSUB3[ncentbins];
TGraphErrors * gN112SUB3[ncentbins];
TGraphErrors * gN112ASUB3[ncentbins];
TGraphErrors * gN112BSUB3[ncentbins];

TGraphErrors * gN1MCm22SUB2_tight2[ncentbins];
TGraphErrors * gN1MCm18SUB2_tight2[ncentbins];
TGraphErrors * gN1MCm14SUB2_tight2[ncentbins];
TGraphErrors * gN1MCp22SUB2_tight2[ncentbins];
TGraphErrors * gN1MCp18SUB2_tight2[ncentbins];
TGraphErrors * gN1MCp14SUB2_tight2[ncentbins];
TGraphErrors * gN1SUB2_tight2[ncentbins];
TGraphErrors * gN1ASUB2_tight2[ncentbins];
TGraphErrors * gN1BSUB2_tight2[ncentbins];
TGraphErrors * gN112SUB2_tight2[ncentbins];
TGraphErrors * gN112ASUB2_tight2[ncentbins];
TGraphErrors * gN112BSUB2_tight2[ncentbins];

TGraphErrors * gN1MCm22SUB3_tight2[ncentbins];
TGraphErrors * gN1MCm18SUB3_tight2[ncentbins];
TGraphErrors * gN1MCm14SUB3_tight2[ncentbins];
TGraphErrors * gN1MCp22SUB3_tight2[ncentbins];
TGraphErrors * gN1MCp18SUB3_tight2[ncentbins];
TGraphErrors * gN1MCp14SUB3_tight2[ncentbins];
TGraphErrors * gN1SUB3_tight2[ncentbins];
TGraphErrors * gN1ASUB3_tight2[ncentbins];
TGraphErrors * gN1BSUB3_tight2[ncentbins];
TGraphErrors * gN112SUB3_tight2[ncentbins];
TGraphErrors * gN112ASUB3_tight2[ncentbins];
TGraphErrors * gN112BSUB3_tight2[ncentbins];

TGraphErrors * gN1MCm22SUB2_narrow[ncentbins];
TGraphErrors * gN1MCm18SUB2_narrow[ncentbins];
TGraphErrors * gN1MCm14SUB2_narrow[ncentbins];
TGraphErrors * gN1MCp22SUB2_narrow[ncentbins];
TGraphErrors * gN1MCp18SUB2_narrow[ncentbins];
TGraphErrors * gN1MCp14SUB2_narrow[ncentbins];
TGraphErrors * gN1SUB2_narrow[ncentbins];
TGraphErrors * gN1ASUB2_narrow[ncentbins];
TGraphErrors * gN1BSUB2_narrow[ncentbins];
TGraphErrors * gN112SUB2_narrow[ncentbins];
TGraphErrors * gN112ASUB2_narrow[ncentbins];
TGraphErrors * gN112BSUB2_narrow[ncentbins];

TGraphErrors * gN1MCm22SUB3_narrow[ncentbins];
TGraphErrors * gN1MCm18SUB3_narrow[ncentbins];
TGraphErrors * gN1MCm14SUB3_narrow[ncentbins];
TGraphErrors * gN1MCp22SUB3_narrow[ncentbins];
TGraphErrors * gN1MCp18SUB3_narrow[ncentbins];
TGraphErrors * gN1MCp14SUB3_narrow[ncentbins];
TGraphErrors * gN1SUB3_narrow[ncentbins];
TGraphErrors * gN1ASUB3_narrow[ncentbins];
TGraphErrors * gN1BSUB3_narrow[ncentbins];
TGraphErrors * gN112SUB3_narrow[ncentbins];
TGraphErrors * gN112ASUB3_narrow[ncentbins];
TGraphErrors * gN112BSUB3_narrow[ncentbins];

TGraphErrors * gN1MCm22SUB2_wide[ncentbins];
TGraphErrors * gN1MCm18SUB2_wide[ncentbins];
TGraphErrors * gN1MCm14SUB2_wide[ncentbins];
TGraphErrors * gN1MCp22SUB2_wide[ncentbins];
TGraphErrors * gN1MCp18SUB2_wide[ncentbins];
TGraphErrors * gN1MCp14SUB2_wide[ncentbins];
TGraphErrors * gN1SUB2_wide[ncentbins];
TGraphErrors * gN1ASUB2_wide[ncentbins];
TGraphErrors * gN1BSUB2_wide[ncentbins];
TGraphErrors * gN112SUB2_wide[ncentbins];
TGraphErrors * gN112ASUB2_wide[ncentbins];
TGraphErrors * gN112BSUB2_wide[ncentbins];

TGraphErrors * gN1MCm22SUB3_wide[ncentbins];
TGraphErrors * gN1MCm18SUB3_wide[ncentbins];
TGraphErrors * gN1MCm14SUB3_wide[ncentbins];
TGraphErrors * gN1MCp22SUB3_wide[ncentbins];
TGraphErrors * gN1MCp18SUB3_wide[ncentbins];
TGraphErrors * gN1MCp14SUB3_wide[ncentbins];
TGraphErrors * gN1SUB3_wide[ncentbins];
TGraphErrors * gN1ASUB3_wide[ncentbins];
TGraphErrors * gN1BSUB3_wide[ncentbins];
TGraphErrors * gN112SUB3_wide[ncentbins];
TGraphErrors * gN112ASUB3_wide[ncentbins];
TGraphErrors * gN112BSUB3_wide[ncentbins];

TH1D * hN1MCm22SUB2[ncentbins];
TH1D * hN1MCm18SUB2[ncentbins];
TH1D * hN1MCm14SUB2[ncentbins];
TH1D * hN1MCp22SUB2[ncentbins];
TH1D * hN1MCp18SUB2[ncentbins];
TH1D * hN1MCp14SUB2[ncentbins];
TH1D * hN1SUB2[ncentbins];
TH1D * hN1ASUB2[ncentbins];
TH1D * hN1BSUB2[ncentbins];
TH1D * hN112SUB2[ncentbins];
TH1D * hN112ASUB2[ncentbins];
TH1D * hN112BSUB2[ncentbins];

TH1D * hN1MCm22SUB3[ncentbins];
TH1D * hN1MCm18SUB3[ncentbins];
TH1D * hN1MCm14SUB3[ncentbins];
TH1D * hN1MCp22SUB3[ncentbins];
TH1D * hN1MCp18SUB3[ncentbins];
TH1D * hN1MCp14SUB3[ncentbins];
TH1D * hN1SUB3[ncentbins];
TH1D * hN1ASUB3[ncentbins];
TH1D * hN1BSUB3[ncentbins];
TH1D * hN112SUB3[ncentbins];
TH1D * hN112ASUB3[ncentbins];
TH1D * hN112BSUB3[ncentbins];

TH1D * hN1MCm22SUB2_tight2[ncentbins];
TH1D * hN1MCm18SUB2_tight2[ncentbins];
TH1D * hN1MCm14SUB2_tight2[ncentbins];
TH1D * hN1MCp22SUB2_tight2[ncentbins];
TH1D * hN1MCp18SUB2_tight2[ncentbins];
TH1D * hN1MCp14SUB2_tight2[ncentbins];
TH1D * hN1SUB2_tight2[ncentbins];
TH1D * hN1ASUB2_tight2[ncentbins];
TH1D * hN1BSUB2_tight2[ncentbins];
TH1D * hN112SUB2_tight2[ncentbins];
TH1D * hN112ASUB2_tight2[ncentbins];
TH1D * hN112BSUB2_tight2[ncentbins];

TH1D * hN1MCm22SUB3_tight2[ncentbins];
TH1D * hN1MCm18SUB3_tight2[ncentbins];
TH1D * hN1MCm14SUB3_tight2[ncentbins];
TH1D * hN1MCp22SUB3_tight2[ncentbins];
TH1D * hN1MCp18SUB3_tight2[ncentbins];
TH1D * hN1MCp14SUB3_tight2[ncentbins];
TH1D * hN1SUB3_tight2[ncentbins];
TH1D * hN1ASUB3_tight2[ncentbins];
TH1D * hN1BSUB3_tight2[ncentbins];
TH1D * hN112SUB3_tight2[ncentbins];
TH1D * hN112ASUB3_tight2[ncentbins];
TH1D * hN112BSUB3_tight2[ncentbins];

TH1D * hN1MCm22SUB2_narrow[ncentbins];
TH1D * hN1MCm18SUB2_narrow[ncentbins];
TH1D * hN1MCm14SUB2_narrow[ncentbins];
TH1D * hN1MCp22SUB2_narrow[ncentbins];
TH1D * hN1MCp18SUB2_narrow[ncentbins];
TH1D * hN1MCp14SUB2_narrow[ncentbins];
TH1D * hN1SUB2_narrow[ncentbins];
TH1D * hN1ASUB2_narrow[ncentbins];
TH1D * hN1BSUB2_narrow[ncentbins];
TH1D * hN112SUB2_narrow[ncentbins];
TH1D * hN112ASUB2_narrow[ncentbins];
TH1D * hN112BSUB2_narrow[ncentbins];

TH1D * hN1MCm22SUB3_narrow[ncentbins];
TH1D * hN1MCm18SUB3_narrow[ncentbins];
TH1D * hN1MCm14SUB3_narrow[ncentbins];
TH1D * hN1MCp22SUB3_narrow[ncentbins];
TH1D * hN1MCp18SUB3_narrow[ncentbins];
TH1D * hN1MCp14SUB3_narrow[ncentbins];
TH1D * hN1SUB3_narrow[ncentbins];
TH1D * hN1ASUB3_narrow[ncentbins];
TH1D * hN1BSUB3_narrow[ncentbins];
TH1D * hN112SUB3_narrow[ncentbins];
TH1D * hN112ASUB3_narrow[ncentbins];
TH1D * hN112BSUB3_narrow[ncentbins];

TH1D * hN1MCm22SUB2_wide[ncentbins];
TH1D * hN1MCm18SUB2_wide[ncentbins];
TH1D * hN1MCm14SUB2_wide[ncentbins];
TH1D * hN1MCp22SUB2_wide[ncentbins];
TH1D * hN1MCp18SUB2_wide[ncentbins];
TH1D * hN1MCp14SUB2_wide[ncentbins];
TH1D * hN1SUB2_wide[ncentbins];
TH1D * hN1ASUB2_wide[ncentbins];
TH1D * hN1BSUB2_wide[ncentbins];
TH1D * hN112SUB2_wide[ncentbins];
TH1D * hN112ASUB2_wide[ncentbins];
TH1D * hN112BSUB2_wide[ncentbins];

TH1D * hN1MCm22SUB3_wide[ncentbins];
TH1D * hN1MCm18SUB3_wide[ncentbins];
TH1D * hN1MCm14SUB3_wide[ncentbins];
TH1D * hN1MCp22SUB3_wide[ncentbins];
TH1D * hN1MCp18SUB3_wide[ncentbins];
TH1D * hN1MCp14SUB3_wide[ncentbins];
TH1D * hN1SUB3_wide[ncentbins];
TH1D * hN1ASUB3_wide[ncentbins];
TH1D * hN1BSUB3_wide[ncentbins];
TH1D * hN112SUB3_wide[ncentbins];
TH1D * hN112ASUB3_wide[ncentbins];
TH1D * hN112BSUB3_wide[ncentbins];

TFile * tfin_nominal;
TFile * tfin_tight2;
TFile * tfin_narrow;
TFile * tfin_wide;

TFile * tfout;

void GetEtaDists()
{
    tfout = new TFile("results/MH_combined_Eta.root","recreate");
    TDirectory * tdir = (TDirectory *) tfout->mkdir("nominal");
    TDirectory * tdir_tight2 = (TDirectory *) tfout->mkdir("tight2");
    TDirectory * tdir_narrow = (TDirectory *) tfout->mkdir("narrow");
    TDirectory * tdir_wide = (TDirectory *) tfout->mkdir("wide");

    for (int i = 0; i<NANALS; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString mtag = Form("%s/%d_%d",ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
            //cout<<mtag.Data()<<endl;
            if (i == 0) gN1MCm22SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 1) gN1MCm18SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 2) gN1MCm14SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 3) gN1MCp22SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 4) gN1MCp18SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 5) gN1MCp14SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 6) gN1MCm22SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 7) gN1MCm18SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 8) gN1MCm14SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 9) gN1MCp22SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 10) gN1MCp18SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 11) gN1MCp14SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 12) gN1SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 13) gN1ASUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 14) gN1BSUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 15) gN112SUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 16) gN112ASUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 17) gN112BSUB2[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 18) gN1SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 19) gN1ASUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 20) gN1BSUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 21) gN112SUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 22) gN112ASUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 23) gN112BSUB3[cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));


            if (i == 0) gN1MCm22SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 1) gN1MCm18SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 2) gN1MCm14SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 3) gN1MCp22SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 4) gN1MCp18SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 5) gN1MCp14SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 6) gN1MCm22SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 7) gN1MCm18SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 8) gN1MCm14SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 9) gN1MCp22SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 10) gN1MCp18SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 11) gN1MCp14SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 12) gN1SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 13) gN1ASUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 14) gN1BSUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 15) gN112SUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 16) gN112ASUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 17) gN112BSUB2_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 18) gN1SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 19) gN1ASUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 20) gN1BSUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 21) gN112SUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 22) gN112ASUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 23) gN112BSUB3_tight2[cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));


            if (i == 0) gN1MCm22SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 1) gN1MCm18SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 2) gN1MCm14SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 3) gN1MCp22SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 4) gN1MCp18SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 5) gN1MCp14SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 6) gN1MCm22SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 7) gN1MCm18SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 8) gN1MCm14SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 9) gN1MCp22SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 10) gN1MCp18SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 11) gN1MCp14SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 12) gN1SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 13) gN1ASUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 14) gN1BSUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 15) gN112SUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 16) gN112ASUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 17) gN112BSUB2_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 18) gN1SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 19) gN1ASUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 20) gN1BSUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 21) gN112SUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 22) gN112ASUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 23) gN112BSUB3_narrow[cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));


            if (i == 0) gN1MCm22SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 1) gN1MCm18SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 2) gN1MCm14SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 3) gN1MCp22SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 4) gN1MCp18SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 5) gN1MCp14SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 6) gN1MCm22SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 7) gN1MCm18SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 8) gN1MCm14SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 9) gN1MCp22SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 10) gN1MCp18SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 11) gN1MCp14SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 12) gN1SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 13) gN1ASUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 14) gN1BSUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 15) gN112SUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 16) gN112ASUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 17) gN112BSUB2_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));

            if (i == 18) gN1SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 19) gN1ASUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 20) gN1BSUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 21) gN112SUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 22) gN112ASUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
            if (i == 23) gN112BSUB3_wide[cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/A_%s",mtag.data(),ANAL[i].Data()));
        }
    }


    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TDirectory * tdcbin = (TDirectory *) tdir->mkdir(Form("%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        tdcbin->cd();

        hN1MCm22SUB2[cbin] = new TH1D(Form("N1MCm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        hN1MCm22SUB2[cbin]->SetStats(0);
        hN1MCm22SUB2[cbin]->SetXTitle("#eta");
        hN1MCm22SUB2[cbin]->SetYTitle("v_{1}");

        hN1MCm18SUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCm14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCp22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCp18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCp14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1MCm18SUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCm14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCp22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCp18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1MCp14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1ASUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1BSUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N112SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N112ASUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB2[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N112BSUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1ASUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N1BSUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N112SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB3[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N112ASUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB23[cbin] = (TH1D *) hN1MCm22SUB2[cbin]->Clone(Form("N112BSUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        GraphToHist( gN1MCm22SUB2[cbin], hN1MCm22SUB2[cbin] );
        GraphToHist( gN1MCm18SUB2[cbin], hN1MCp18SUB2[cbin] );
        GraphToHist( gN1MCm14SUB2[cbin], hN1MCp14SUB2[cbin] );
        GraphToHist( gN1MCp22SUB2[cbin], hN1MCp22SUB2[cbin] );
        GraphToHist( gN1MCp18SUB2[cbin], hN1MCp18SUB2[cbin] );
        GraphToHist( gN1MCp14SUB2[cbin], hN1MCp14SUB2[cbin] );

        GraphToHist( gN1MCm22SUB3[cbin], hN1MCm22SUB3[cbin] );
        GraphToHist( gN1MCm18SUB3[cbin], hN1MCp18SUB3[cbin] );
        GraphToHist( gN1MCm14SUB3[cbin], hN1MCp14SUB3[cbin] );
        GraphToHist( gN1MCp22SUB3[cbin], hN1MCp22SUB3[cbin] );
        GraphToHist( gN1MCp18SUB3[cbin], hN1MCp18SUB3[cbin] );
        GraphToHist( gN1MCp14SUB3[cbin], hN1MCp14SUB3[cbin] );

        GraphToHist( gN1SUB2[cbin], hN1SUB2[cbin] );
        GraphToHist( gN1ASUB2[cbin], hN1ASUB2[cbin] );
        GraphToHist( gN1BSUB2[cbin], hN1BSUB2[cbin] );
        GraphToHist( gN112SUB2[cbin], hN112SUB2[cbin] );
        GraphToHist( gN112ASUB2[cbin], hN112ASUB2[cbin] );
        GraphToHist( gN112BSUB2[cbin], hN112BSUB2[cbin] );

        GraphToHist( gN1SUB3[cbin], hN1SUB3[cbin] );
        GraphToHist( gN1ASUB3[cbin], hN1ASUB3[cbin] );
        GraphToHist( gN1BSUB3[cbin], hN1BSUB3[cbin] );
        GraphToHist( gN112SUB3[cbin], hN112SUB3[cbin] );
        GraphToHist( gN112ASUB3[cbin], hN112ASUB3[cbin] );
        GraphToHist( gN112BSUB3[cbin], hN112BSUB3[cbin] );

        hN1MCm22SUB2[cbin]->Write();
        hN1MCm18SUB2[cbin]->Write();
        hN1MCm14SUB2[cbin]->Write();
        hN1MCp22SUB2[cbin]->Write();
        hN1MCp18SUB2[cbin]->Write();
        hN1MCp14SUB2[cbin]->Write();

        hN1MCm22SUB3[cbin]->Write();
        hN1MCm18SUB3[cbin]->Write();
        hN1MCm14SUB3[cbin]->Write();
        hN1MCp22SUB3[cbin]->Write();
        hN1MCp18SUB3[cbin]->Write();
        hN1MCp14SUB3[cbin]->Write();

        hN1SUB2[cbin]->Write();
        hN1ASUB2[cbin]->Write();
        hN1BSUB2[cbin]->Write();
        hN112SUB2[cbin]->Write();
        hN112ASUB2[cbin]->Write();
        hN112BSUB2[cbin]->Write();

        hN1SUB3[cbin]->Write();
        hN1ASUB3[cbin]->Write();
        hN1BSUB3[cbin]->Write();
        hN112SUB3[cbin]->Write();
        hN112ASUB3[cbin]->Write();
        hN112BSUB23[cbin]->Write();
    }

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TDirectory * tdcbin = (TDirectory *) tdir_tight2->mkdir(Form("%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        tdcbin->cd();

        hN1MCm22SUB2_tight2[cbin] = new TH1D(Form("N1MCm22SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        hN1MCm22SUB2_tight2[cbin]->SetStats(0);
        hN1MCm22SUB2_tight2[cbin]->SetXTitle("#eta");
        hN1MCm22SUB2_tight2[cbin]->SetYTitle("v_{1} (tight2 cuts)");

        hN1MCm18SUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCm18SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCm14SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCp22SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCp18SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCp14SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1MCm18SUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCm18SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCm14SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCp22SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCp18SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1MCp14SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1ASUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1BSUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N112SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N112ASUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB2_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N112BSUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1ASUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N1BSUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N112SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N112ASUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB3_tight2[cbin] = (TH1D *) hN1MCm22SUB2_tight2[cbin]->Clone(Form("N112BSUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        GraphToHist( gN1MCm22SUB2_tight2[cbin], hN1MCm22SUB2_tight2[cbin] );
        GraphToHist( gN1MCm18SUB2_tight2[cbin], hN1MCp18SUB2_tight2[cbin] );
        GraphToHist( gN1MCm14SUB2_tight2[cbin], hN1MCp14SUB2_tight2[cbin] );
        GraphToHist( gN1MCp22SUB2_tight2[cbin], hN1MCp22SUB2_tight2[cbin] );
        GraphToHist( gN1MCp18SUB2_tight2[cbin], hN1MCp18SUB2_tight2[cbin] );
        GraphToHist( gN1MCp14SUB2_tight2[cbin], hN1MCp14SUB2_tight2[cbin] );

        GraphToHist( gN1MCm22SUB3_tight2[cbin], hN1MCm22SUB3_tight2[cbin] );
        GraphToHist( gN1MCm18SUB3_tight2[cbin], hN1MCp18SUB3_tight2[cbin] );
        GraphToHist( gN1MCm14SUB3_tight2[cbin], hN1MCp14SUB3_tight2[cbin] );
        GraphToHist( gN1MCp22SUB3_tight2[cbin], hN1MCp22SUB3_tight2[cbin] );
        GraphToHist( gN1MCp18SUB3_tight2[cbin], hN1MCp18SUB3_tight2[cbin] );
        GraphToHist( gN1MCp14SUB3_tight2[cbin], hN1MCp14SUB3_tight2[cbin] );

        GraphToHist( gN1SUB2_tight2[cbin], hN1SUB2_tight2[cbin] );
        GraphToHist( gN1ASUB2_tight2[cbin], hN1ASUB2_tight2[cbin] );
        GraphToHist( gN1BSUB2_tight2[cbin], hN1BSUB2_tight2[cbin] );
        GraphToHist( gN112SUB2_tight2[cbin], hN112SUB2_tight2[cbin] );
        GraphToHist( gN112ASUB2_tight2[cbin], hN112ASUB2_tight2[cbin] );
        GraphToHist( gN112BSUB2_tight2[cbin], hN112BSUB2_tight2[cbin] );

        GraphToHist( gN1SUB3_tight2[cbin], hN1SUB3_tight2[cbin] );
        GraphToHist( gN1ASUB3_tight2[cbin], hN1ASUB3_tight2[cbin] );
        GraphToHist( gN1BSUB3_tight2[cbin], hN1BSUB3_tight2[cbin] );
        GraphToHist( gN112SUB3_tight2[cbin], hN112SUB3_tight2[cbin] );
        GraphToHist( gN112ASUB3_tight2[cbin], hN112ASUB3_tight2[cbin] );
        GraphToHist( gN112BSUB3_tight2[cbin], hN112BSUB3_tight2[cbin] );

        hN1MCm22SUB2_tight2[cbin]->Write();
        hN1MCm18SUB2_tight2[cbin]->Write();
        hN1MCm14SUB2_tight2[cbin]->Write();
        hN1MCp22SUB2_tight2[cbin]->Write();
        hN1MCp18SUB2_tight2[cbin]->Write();
        hN1MCp14SUB2_tight2[cbin]->Write();

        hN1MCm22SUB3_tight2[cbin]->Write();
        hN1MCm18SUB3_tight2[cbin]->Write();
        hN1MCm14SUB3_tight2[cbin]->Write();
        hN1MCp22SUB3_tight2[cbin]->Write();
        hN1MCp18SUB3_tight2[cbin]->Write();
        hN1MCp14SUB3_tight2[cbin]->Write();

        hN1SUB2_tight2[cbin]->Write();
        hN1ASUB2_tight2[cbin]->Write();
        hN1BSUB2_tight2[cbin]->Write();
        hN112SUB2_tight2[cbin]->Write();
        hN112ASUB2_tight2[cbin]->Write();
        hN112BSUB2_tight2[cbin]->Write();

        hN1SUB3_tight2[cbin]->Write();
        hN1ASUB3_tight2[cbin]->Write();
        hN1BSUB3_tight2[cbin]->Write();
        hN112SUB3_tight2[cbin]->Write();
        hN112ASUB3_tight2[cbin]->Write();
        hN112BSUB3_tight2[cbin]->Write();
    }

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TDirectory * tdcbin = (TDirectory *) tdir_narrow->mkdir(Form("%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        tdcbin->cd();

        hN1MCm22SUB2_narrow[cbin] = new TH1D(Form("N1MCm22SUB2_narrow_%s_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
        hN1MCm22SUB2_narrow[cbin]->SetStats(0);
        hN1MCm22SUB2_narrow[cbin]->SetXTitle("#eta");
        hN1MCm22SUB2_narrow[cbin]->SetYTitle("v_{1} (narrow cuts)");

        hN1MCm18SUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCm18SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCm14SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCp22SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCp18SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCp14SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1MCm18SUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCm18SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCm14SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCp22SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCp18SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1MCp14SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1ASUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1BSUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N112SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N112ASUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB2_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N112BSUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1ASUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N1BSUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N112SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N112ASUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB3_narrow[cbin] = (TH1D *) hN1MCm22SUB2_narrow[cbin]->Clone(Form("N112BSUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        GraphToHist( gN1MCm22SUB2_narrow[cbin], hN1MCm22SUB2_narrow[cbin] );
        GraphToHist( gN1MCm18SUB2_narrow[cbin], hN1MCp18SUB2_narrow[cbin] );
        GraphToHist( gN1MCm14SUB2_narrow[cbin], hN1MCp14SUB2_narrow[cbin] );
        GraphToHist( gN1MCp22SUB2_narrow[cbin], hN1MCp22SUB2_narrow[cbin] );
        GraphToHist( gN1MCp18SUB2_narrow[cbin], hN1MCp18SUB2_narrow[cbin] );
        GraphToHist( gN1MCp14SUB2_narrow[cbin], hN1MCp14SUB2_narrow[cbin] );

        GraphToHist( gN1MCm22SUB3_narrow[cbin], hN1MCm22SUB3_narrow[cbin] );
        GraphToHist( gN1MCm18SUB3_narrow[cbin], hN1MCp18SUB3_narrow[cbin] );
        GraphToHist( gN1MCm14SUB3_narrow[cbin], hN1MCp14SUB3_narrow[cbin] );
        GraphToHist( gN1MCp22SUB3_narrow[cbin], hN1MCp22SUB3_narrow[cbin] );
        GraphToHist( gN1MCp18SUB3_narrow[cbin], hN1MCp18SUB3_narrow[cbin] );
        GraphToHist( gN1MCp14SUB3_narrow[cbin], hN1MCp14SUB3_narrow[cbin] );

        GraphToHist( gN1SUB2_narrow[cbin], hN1SUB2_narrow[cbin] );
        GraphToHist( gN1ASUB2_narrow[cbin], hN1ASUB2_narrow[cbin] );
        GraphToHist( gN1BSUB2_narrow[cbin], hN1BSUB2_narrow[cbin] );
        GraphToHist( gN112SUB2_narrow[cbin], hN112SUB2_narrow[cbin] );
        GraphToHist( gN112ASUB2_narrow[cbin], hN112ASUB2_narrow[cbin] );
        GraphToHist( gN112BSUB2_narrow[cbin], hN112BSUB2_narrow[cbin] );

        GraphToHist( gN1SUB3_narrow[cbin], hN1SUB3_narrow[cbin] );
        GraphToHist( gN1ASUB3_narrow[cbin], hN1ASUB3_narrow[cbin] );
        GraphToHist( gN1BSUB3_narrow[cbin], hN1BSUB3_narrow[cbin] );
        GraphToHist( gN112SUB3_narrow[cbin], hN112SUB3_narrow[cbin] );
        GraphToHist( gN112ASUB3_narrow[cbin], hN112ASUB3_narrow[cbin] );
        GraphToHist( gN112BSUB3_narrow[cbin], hN112BSUB3_narrow[cbin] );

        hN1MCm22SUB2_narrow[cbin]->Write();
        hN1MCm18SUB2_narrow[cbin]->Write();
        hN1MCm14SUB2_narrow[cbin]->Write();
        hN1MCp22SUB2_narrow[cbin]->Write();
        hN1MCp18SUB2_narrow[cbin]->Write();
        hN1MCp14SUB2_narrow[cbin]->Write();

        hN1MCm22SUB3_narrow[cbin]->Write();
        hN1MCm18SUB3_narrow[cbin]->Write();
        hN1MCm14SUB3_narrow[cbin]->Write();
        hN1MCp22SUB3_narrow[cbin]->Write();
        hN1MCp18SUB3_narrow[cbin]->Write();
        hN1MCp14SUB3_narrow[cbin]->Write();

        hN1SUB2_narrow[cbin]->Write();
        hN1ASUB2_narrow[cbin]->Write();
        hN1BSUB2_narrow[cbin]->Write();
        hN112SUB2_narrow[cbin]->Write();
        hN112ASUB2_narrow[cbin]->Write();
        hN112BSUB2_narrow[cbin]->Write();

        hN1SUB3_narrow[cbin]->Write();
        hN1ASUB3_narrow[cbin]->Write();
        hN1BSUB3_narrow[cbin]->Write();
        hN112SUB3_narrow[cbin]->Write();
        hN112ASUB3_narrow[cbin]->Write();
        hN112BSUB3_narrow[cbin]->Write();
    }

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TDirectory * tdcbin = (TDirectory *) tdir_wide->mkdir(Form("%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        tdcbin->cd();

        hN1MCm22SUB2_wide[cbin] = new TH1D(Form("N1MCm22SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
        hN1MCm22SUB2_wide[cbin]->SetStats(0);
        hN1MCm22SUB2_wide[cbin]->SetXTitle("#eta");
        hN1MCm22SUB2_wide[cbin]->SetYTitle("v_{1} (wide cuts)");

        hN1MCm18SUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCm18SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCm14SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCp22SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCp18SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCp14SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1MCm18SUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCm18SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCm14SUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCm14SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp22SUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCp22SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp18SUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCp18SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCp14SUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1MCp14SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1ASUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1BSUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N112SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N112ASUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB2_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N112BSUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        hN1SUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1ASUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1ASUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1BSUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N1BSUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112SUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N112SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112ASUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N112ASUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN112BSUB3_wide[cbin] = (TH1D *) hN1MCm22SUB2_wide[cbin]->Clone(Form("N112BSUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));

        GraphToHist( gN1MCm22SUB2_wide[cbin], hN1MCm22SUB2_wide[cbin] );
        GraphToHist( gN1MCm18SUB2_wide[cbin], hN1MCp18SUB2_wide[cbin] );
        GraphToHist( gN1MCm14SUB2_wide[cbin], hN1MCp14SUB2_wide[cbin] );
        GraphToHist( gN1MCp22SUB2_wide[cbin], hN1MCp22SUB2_wide[cbin] );
        GraphToHist( gN1MCp18SUB2_wide[cbin], hN1MCp18SUB2_wide[cbin] );
        GraphToHist( gN1MCp14SUB2_wide[cbin], hN1MCp14SUB2_wide[cbin] );

        GraphToHist( gN1MCm22SUB3_wide[cbin], hN1MCm22SUB3_wide[cbin] );
        GraphToHist( gN1MCm18SUB3_wide[cbin], hN1MCp18SUB3_wide[cbin] );
        GraphToHist( gN1MCm14SUB3_wide[cbin], hN1MCp14SUB3_wide[cbin] );
        GraphToHist( gN1MCp22SUB3_wide[cbin], hN1MCp22SUB3_wide[cbin] );
        GraphToHist( gN1MCp18SUB3_wide[cbin], hN1MCp18SUB3_wide[cbin] );
        GraphToHist( gN1MCp14SUB3_wide[cbin], hN1MCp14SUB3_wide[cbin] );

        GraphToHist( gN1SUB2_wide[cbin], hN1SUB2_wide[cbin] );
        GraphToHist( gN1ASUB2_wide[cbin], hN1ASUB2_wide[cbin] );
        GraphToHist( gN1BSUB2_wide[cbin], hN1BSUB2_wide[cbin] );
        GraphToHist( gN112SUB2_wide[cbin], hN112SUB2_wide[cbin] );
        GraphToHist( gN112ASUB2_wide[cbin], hN112ASUB2_wide[cbin] );
        GraphToHist( gN112BSUB2_wide[cbin], hN112BSUB2_wide[cbin] );

        GraphToHist( gN1SUB3_wide[cbin], hN1SUB3_wide[cbin] );
        GraphToHist( gN1ASUB3_wide[cbin], hN1ASUB3_wide[cbin] );
        GraphToHist( gN1BSUB3_wide[cbin], hN1BSUB3_wide[cbin] );
        GraphToHist( gN112SUB3_wide[cbin], hN112SUB3_wide[cbin] );
        GraphToHist( gN112ASUB3_wide[cbin], hN112ASUB3_wide[cbin] );
        GraphToHist( gN112BSUB3_wide[cbin], hN112BSUB3_wide[cbin] );

        hN1MCm22SUB2_wide[cbin]->Write();
        hN1MCm18SUB2_wide[cbin]->Write();
        hN1MCm14SUB2_wide[cbin]->Write();
        hN1MCp22SUB2_wide[cbin]->Write();
        hN1MCp18SUB2_wide[cbin]->Write();
        hN1MCp14SUB2_wide[cbin]->Write();

        hN1MCm22SUB3_wide[cbin]->Write();
        hN1MCm18SUB3_wide[cbin]->Write();
        hN1MCm14SUB3_wide[cbin]->Write();
        hN1MCp22SUB3_wide[cbin]->Write();
        hN1MCp18SUB3_wide[cbin]->Write();
        hN1MCp14SUB3_wide[cbin]->Write();

        hN1SUB2_wide[cbin]->Write();
        hN1ASUB2_wide[cbin]->Write();
        hN1BSUB2_wide[cbin]->Write();
        hN112SUB2_wide[cbin]->Write();
        hN112ASUB2_wide[cbin]->Write();
        hN112BSUB2_wide[cbin]->Write();

        hN1SUB3_wide[cbin]->Write();
        hN1ASUB3_wide[cbin]->Write();
        hN1BSUB3_wide[cbin]->Write();
        hN112SUB3_wide[cbin]->Write();
        hN112ASUB3_wide[cbin]->Write();
        hN112BSUB3_wide[cbin]->Write();
    }

    tfout->Close();

}

void combineEta() {

    TH1::SetDefaultSumw2();

    if (!fopen("results","r")) system("mkdir results");
    tfout = new TFile("results/MH_combined_Pt.root","recreate");

    GetPtDists( "MH_hists" );
    GetPtDists( "MH_tight2_hists" );
    GetPtDists( "MH_narrow_hists" );
    GetPtDists( "MH_wide_hists" );

    cout << "Pt distributions saved to results/MH_combined_Pt.root" << endl;
    tfout->Close();

}
