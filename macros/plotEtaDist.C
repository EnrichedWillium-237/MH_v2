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
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCp22SUB2", "N1MCp18SUB2", "N1MCp14SUB2",
                  "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCp22SUB3", "N1MCp18SUB3", "N1MCp14SUB3",
                  "N1SUB2",      "N1ASUB2",     "N1BSUB2",     "N112SUB2",    "N112ASUB2",   "N112BSUB2",
                  "N1SUB3",      "N1ASUB3",     "N1BSUB3",     "N112SUB3",    "N112ASUB3",   "N112BSUB3"
};

TH1D * N1MCm22SUB2[ncentbins];
TH1D * N1MCm18SUB2[ncentbins];
TH1D * N1MCm14SUB2[ncentbins];
TH1D * N1MCp22SUB2[ncentbins];
TH1D * N1MCp18SUB2[ncentbins];
TH1D * N1MCp14SUB2[ncentbins];

TH1D * N1MCm22SUB3[ncentbins];
TH1D * N1MCm18SUB3[ncentbins];
TH1D * N1MCm14SUB3[ncentbins];
TH1D * N1MCp22SUB3[ncentbins];
TH1D * N1MCp18SUB3[ncentbins];
TH1D * N1MCp14SUB3[ncentbins];

TH1D * N1SUB2[ncentbins];
TH1D * N1ASUB2[ncentbins];
TH1D * N1BSUB2[ncentbins];
TH1D * N112SUB2[ncentbins];
TH1D * N112ASUB2[ncentbins];
TH1D * N112BSUB2[ncentbins];

TH1D * N1SUB3[ncentbins];
TH1D * N1ASUB3[ncentbins];
TH1D * N1BSUB3[ncentbins];
TH1D * N112SUB3[ncentbins];
TH1D * N112ASUB3[ncentbins];
TH1D * N112BSUB3[ncentbins];


TH1D * N1MCm22SUB2_tight2[ncentbins];
TH1D * N1MCm18SUB2_tight2[ncentbins];
TH1D * N1MCm14SUB2_tight2[ncentbins];
TH1D * N1MCp22SUB2_tight2[ncentbins];
TH1D * N1MCp18SUB2_tight2[ncentbins];
TH1D * N1MCp14SUB2_tight2[ncentbins];

TH1D * N1MCm22SUB3_tight2[ncentbins];
TH1D * N1MCm18SUB3_tight2[ncentbins];
TH1D * N1MCm14SUB3_tight2[ncentbins];
TH1D * N1MCp22SUB3_tight2[ncentbins];
TH1D * N1MCp18SUB3_tight2[ncentbins];
TH1D * N1MCp14SUB3_tight2[ncentbins];

TH1D * N1SUB2_tight2[ncentbins];
TH1D * N1ASUB2_tight2[ncentbins];
TH1D * N1BSUB2_tight2[ncentbins];
TH1D * N112SUB2_tight2[ncentbins];
TH1D * N112ASUB2_tight2[ncentbins];
TH1D * N112BSUB2_tight2[ncentbins];

TH1D * N1SUB3_tight2[ncentbins];
TH1D * N1ASUB3_tight2[ncentbins];
TH1D * N1BSUB3_tight2[ncentbins];
TH1D * N112SUB3_tight2[ncentbins];
TH1D * N112ASUB3_tight2[ncentbins];
TH1D * N112BSUB3_tight2[ncentbins];


TH1D * N1MCm22SUB2_wide[ncentbins];
TH1D * N1MCm18SUB2_wide[ncentbins];
TH1D * N1MCm14SUB2_wide[ncentbins];
TH1D * N1MCp22SUB2_wide[ncentbins];
TH1D * N1MCp18SUB2_wide[ncentbins];
TH1D * N1MCp14SUB2_wide[ncentbins];

TH1D * N1MCm22SUB3_wide[ncentbins];
TH1D * N1MCm18SUB3_wide[ncentbins];
TH1D * N1MCm14SUB3_wide[ncentbins];
TH1D * N1MCp22SUB3_wide[ncentbins];
TH1D * N1MCp18SUB3_wide[ncentbins];
TH1D * N1MCp14SUB3_wide[ncentbins];

TH1D * N1SUB2_wide[ncentbins];
TH1D * N1ASUB2_wide[ncentbins];
TH1D * N1BSUB2_wide[ncentbins];
TH1D * N112SUB2_wide[ncentbins];
TH1D * N112ASUB2_wide[ncentbins];
TH1D * N112BSUB2_wide[ncentbins];

TH1D * N1SUB3_wide[ncentbins];
TH1D * N1ASUB3_wide[ncentbins];
TH1D * N1BSUB3_wide[ncentbins];
TH1D * N112SUB3_wide[ncentbins];
TH1D * N112ASUB3_wide[ncentbins];
TH1D * N112BSUB3_wide[ncentbins];


TH1D * N1MCm22SUB2_narrow[ncentbins];
TH1D * N1MCm18SUB2_narrow[ncentbins];
TH1D * N1MCm14SUB2_narrow[ncentbins];
TH1D * N1MCp22SUB2_narrow[ncentbins];
TH1D * N1MCp18SUB2_narrow[ncentbins];
TH1D * N1MCp14SUB2_narrow[ncentbins];

TH1D * N1MCm22SUB3_narrow[ncentbins];
TH1D * N1MCm18SUB3_narrow[ncentbins];
TH1D * N1MCm14SUB3_narrow[ncentbins];
TH1D * N1MCp22SUB3_narrow[ncentbins];
TH1D * N1MCp18SUB3_narrow[ncentbins];
TH1D * N1MCp14SUB3_narrow[ncentbins];

TH1D * N1SUB2_narrow[ncentbins];
TH1D * N1ASUB2_narrow[ncentbins];
TH1D * N1BSUB2_narrow[ncentbins];
TH1D * N112SUB2_narrow[ncentbins];
TH1D * N112ASUB2_narrow[ncentbins];
TH1D * N112BSUB2_narrow[ncentbins];

TH1D * N1SUB3_narrow[ncentbins];
TH1D * N1ASUB3_narrow[ncentbins];
TH1D * N1BSUB3_narrow[ncentbins];
TH1D * N112SUB3_narrow[ncentbins];
TH1D * N112ASUB3_narrow[ncentbins];
TH1D * N112BSUB3_narrow[ncentbins];

TFile * tfin;

void plotEtaDist()
{
    tfin = new TFile("results/MH_combined_Eta.root","read");

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TString mtag = Form("MH_nominal/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1MCm22SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB3[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));


        mtag = Form("MH_tight2/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1MCm22SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB3_tight2[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));


        mtag = Form("MH_wide/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1MCm22SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB3_wide[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));


        mtag = Form("MH_narrow/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1MCm22SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm18SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm14SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp18SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp14SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp14SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        N1SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB3_narrow[cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
    }

    if (!fopen("results/Eta_Distributions","r")) system("mkdir results/Eta_Distributions");

    TH1D * h0 = new TH1D("h0", "", 100, 0.0, 12.0);
    h0->SetStats(0);
    h0->SetXTitle("p_{T} (GeV/c)");
    h0->SetYTitle("v_{1}");
    h0->GetYaxis()->SetRangeUser(-0.03, 0.2);

    for (int cbin = )

}
