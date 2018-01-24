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

TH1D * N1MC22SUB2[ncentbins];
TH1D * N1MC18SUB2[ncentbins];
TH1D * N1MC14SUB2[ncentbins];
TH1D * N1MC22SUB3[ncentbins];
TH1D * N1MC18SUB3[ncentbins];
TH1D * N1MC14SUB3[ncentbins];


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

TH1D * N1MC22SUB2_tight2[ncentbins];
TH1D * N1MC18SUB2_tight2[ncentbins];
TH1D * N1MC14SUB2_tight2[ncentbins];
TH1D * N1MC22SUB3_tight2[ncentbins];
TH1D * N1MC18SUB3_tight2[ncentbins];
TH1D * N1MC14SUB3_tight2[ncentbins];


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

TH1D * N1MC22SUB2_wide[ncentbins];
TH1D * N1MC18SUB2_wide[ncentbins];
TH1D * N1MC14SUB2_wide[ncentbins];
TH1D * N1MC22SUB3_wide[ncentbins];
TH1D * N1MC18SUB3_wide[ncentbins];
TH1D * N1MC14SUB3_wide[ncentbins];


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

TH1D * N1MC22SUB2_narrow[ncentbins];
TH1D * N1MC18SUB2_narrow[ncentbins];
TH1D * N1MC14SUB2_narrow[ncentbins];
TH1D * N1MC22SUB3_narrow[ncentbins];
TH1D * N1MC18SUB3_narrow[ncentbins];
TH1D * N1MC14SUB3_narrow[ncentbins];

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

        N1MC22SUB2[cbin] = new TH1D(Form("N1MC22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB2[cbin] = new TH1D(Form("N1MC18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB2[cbin] = new TH1D(Form("N1MC14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC22SUB3[cbin] = new TH1D(Form("N1MC22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB3[cbin] = new TH1D(Form("N1MC18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB3[cbin] = new TH1D(Form("N1MC14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);


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

        N1MC22SUB2_tight2[cbin] = new TH1D(Form("N1MC22SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB2_tight2[cbin] = new TH1D(Form("N1MC18SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB2_tight2[cbin] = new TH1D(Form("N1MC14SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC22SUB3_tight2[cbin] = new TH1D(Form("N1MC22SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB3_tight2[cbin] = new TH1D(Form("N1MC18SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB3_tight2[cbin] = new TH1D(Form("N1MC14SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);


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

        N1MC22SUB2_wide[cbin] = new TH1D(Form("N1MC22SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB2_wide[cbin] = new TH1D(Form("N1MC18SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB2_wide[cbin] = new TH1D(Form("N1MC14SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC22SUB3_wide[cbin] = new TH1D(Form("N1MC22SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB3_wide[cbin] = new TH1D(Form("N1MC18SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB3_wide[cbin] = new TH1D(Form("N1MC14SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);


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

        N1MC22SUB2_narrow[cbin] = new TH1D(Form("N1MC22SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB2_narrow[cbin] = new TH1D(Form("N1MC18SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB2_narrow[cbin] = new TH1D(Form("N1MC14SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC22SUB3_narrow[cbin] = new TH1D(Form("N1MC22SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC18SUB3_narrow[cbin] = new TH1D(Form("N1MC18SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        N1MC14SUB3_narrow[cbin] = new TH1D(Form("N1MC14SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
    }

    double emin1 = -2.4;
    double emax1 = -0.001;
    double emin2 = 0.001;
    double emax2 = 2.4;
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        for (int ebin = 1; ebin<=netabins; ebin++) {
            if (N1MC22SUB2[cbin]->GetBinCenter(ebin)>emin1 && N1MC22SUB2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC22SUB2[cbin]->SetBinContent(ebin, N1MCp22SUB2[cbin]->GetBinContent(ebin));
                N1MC22SUB2[cbin]->SetBinError(ebin, N1MCp22SUB2[cbin]->GetBinError(ebin));
            } else {
                N1MC22SUB2[cbin]->SetBinContent(ebin, N1MCm22SUB2[cbin]->GetBinContent(ebin));
                N1MC22SUB2[cbin]->SetBinError(ebin, N1MCm22SUB2[cbin]->GetBinError(ebin));
            }
            if (N1MC18SUB2[cbin]->GetBinCenter(ebin)>emin1 && N1MC18SUB2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC18SUB2[cbin]->SetBinContent(ebin, N1MCp18SUB2[cbin]->GetBinContent(ebin));
                N1MC18SUB2[cbin]->SetBinError(ebin, N1MCp18SUB2[cbin]->GetBinError(ebin));
            } else {
                N1MC18SUB2[cbin]->SetBinContent(ebin, N1MCm18SUB2[cbin]->GetBinContent(ebin));
                N1MC18SUB2[cbin]->SetBinError(ebin, N1MCm18SUB2[cbin]->GetBinError(ebin));
            }
            if (N1MC14SUB2[cbin]->GetBinCenter(ebin)>emin1 && N1MC14SUB2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC14SUB2[cbin]->SetBinContent(ebin, N1MCp14SUB2[cbin]->GetBinContent(ebin));
                N1MC14SUB2[cbin]->SetBinError(ebin, N1MCp14SUB2[cbin]->GetBinError(ebin));
            } else {
                N1MC14SUB2[cbin]->SetBinContent(ebin, N1MCm14SUB2[cbin]->GetBinContent(ebin));
                N1MC14SUB2[cbin]->SetBinError(ebin, N1MCm14SUB2[cbin]->GetBinError(ebin));
            }
            if (N1MC22SUB3[cbin]->GetBinCenter(ebin)>emin1 && N1MC22SUB3[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC22SUB3[cbin]->SetBinContent(ebin, N1MCp22SUB3[cbin]->GetBinContent(ebin));
                N1MC22SUB3[cbin]->SetBinError(ebin, N1MCp22SUB3[cbin]->GetBinError(ebin));
            } else {
                N1MC22SUB3[cbin]->SetBinContent(ebin, N1MCm22SUB3[cbin]->GetBinContent(ebin));
                N1MC22SUB3[cbin]->SetBinError(ebin, N1MCm22SUB3[cbin]->GetBinError(ebin));
            }
            if (N1MC18SUB3[cbin]->GetBinCenter(ebin)>emin1 && N1MC18SUB3[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC18SUB3[cbin]->SetBinContent(ebin, N1MCp18SUB3[cbin]->GetBinContent(ebin));
                N1MC18SUB3[cbin]->SetBinError(ebin, N1MCp18SUB3[cbin]->GetBinError(ebin));
            } else {
                N1MC18SUB3[cbin]->SetBinContent(ebin, N1MCm18SUB3[cbin]->GetBinContent(ebin));
                N1MC18SUB3[cbin]->SetBinError(ebin, N1MCm18SUB3[cbin]->GetBinError(ebin));
            }
            if (N1MC14SUB3[cbin]->GetBinCenter(ebin)>emin1 && N1MC14SUB3[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC14SUB3[cbin]->SetBinContent(ebin, N1MCp14SUB3[cbin]->GetBinContent(ebin));
                N1MC14SUB3[cbin]->SetBinError(ebin, N1MCp14SUB3[cbin]->GetBinError(ebin));
            } else {
                N1MC14SUB3[cbin]->SetBinContent(ebin, N1MCm14SUB3[cbin]->GetBinContent(ebin));
                N1MC14SUB3[cbin]->SetBinError(ebin, N1MCm14SUB3[cbin]->GetBinError(ebin));
            }

            if (N1MC22SUB2_tight2[cbin]->GetBinCenter(ebin)>emin1 && N1MC22SUB2_tight2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC22SUB2_tight2[cbin]->SetBinContent(ebin, N1MCp22SUB2_tight2[cbin]->GetBinContent(ebin));
                N1MC22SUB2_tight2[cbin]->SetBinError(ebin, N1MCp22SUB2_tight2[cbin]->GetBinError(ebin));
            } else {
                N1MC22SUB2_tight2[cbin]->SetBinContent(ebin, N1MCm22SUB2_tight2[cbin]->GetBinContent(ebin));
                N1MC22SUB2_tight2[cbin]->SetBinError(ebin, N1MCm22SUB2_tight2[cbin]->GetBinError(ebin));
            }
            if (N1MC18SUB2_tight2[cbin]->GetBinCenter(ebin)>emin1 && N1MC18SUB2_tight2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC18SUB2_tight2[cbin]->SetBinContent(ebin, N1MCp18SUB2_tight2[cbin]->GetBinContent(ebin));
                N1MC18SUB2_tight2[cbin]->SetBinError(ebin, N1MCp18SUB2_tight2[cbin]->GetBinError(ebin));
            } else {
                N1MC18SUB2_tight2[cbin]->SetBinContent(ebin, N1MCm18SUB2_tight2[cbin]->GetBinContent(ebin));
                N1MC18SUB2_tight2[cbin]->SetBinError(ebin, N1MCm18SUB2_tight2[cbin]->GetBinError(ebin));
            }
            if (N1MC14SUB2_tight2[cbin]->GetBinCenter(ebin)>emin1 && N1MC14SUB2_tight2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC14SUB2_tight2[cbin]->SetBinContent(ebin, N1MCp14SUB2_tight2[cbin]->GetBinContent(ebin));
                N1MC14SUB2_tight2[cbin]->SetBinError(ebin, N1MCp14SUB2_tight2[cbin]->GetBinError(ebin));
            } else {
                N1MC14SUB2_tight2[cbin]->SetBinContent(ebin, N1MCm14SUB2_tight2[cbin]->GetBinContent(ebin));
                N1MC14SUB2_tight2[cbin]->SetBinError(ebin, N1MCm14SUB2_tight2[cbin]->GetBinError(ebin));
            }
            if (N1MC22SUB3_tight2[cbin]->GetBinCenter(ebin)>emin1 && N1MC22SUB3_tight2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC22SUB3_tight2[cbin]->SetBinContent(ebin, N1MCp22SUB3_tight2[cbin]->GetBinContent(ebin));
                N1MC22SUB3_tight2[cbin]->SetBinError(ebin, N1MCp22SUB3_tight2[cbin]->GetBinError(ebin));
            } else {
                N1MC22SUB3_tight2[cbin]->SetBinContent(ebin, N1MCm22SUB3_tight2[cbin]->GetBinContent(ebin));
                N1MC22SUB3_tight2[cbin]->SetBinError(ebin, N1MCm22SUB3_tight2[cbin]->GetBinError(ebin));
            }
            if (N1MC18SUB3_tight2[cbin]->GetBinCenter(ebin)>emin1 && N1MC18SUB3_tight2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC18SUB3_tight2[cbin]->SetBinContent(ebin, N1MCp18SUB3_tight2[cbin]->GetBinContent(ebin));
                N1MC18SUB3_tight2[cbin]->SetBinError(ebin, N1MCp18SUB3_tight2[cbin]->GetBinError(ebin));
            } else {
                N1MC18SUB3_tight2[cbin]->SetBinContent(ebin, N1MCm18SUB3_tight2[cbin]->GetBinContent(ebin));
                N1MC18SUB3_tight2[cbin]->SetBinError(ebin, N1MCm18SUB3_tight2[cbin]->GetBinError(ebin));
            }
            if (N1MC14SUB3_tight2[cbin]->GetBinCenter(ebin)>emin1 && N1MC14SUB3_tight2[cbin]->GetBinCenter(ebin)<emax1) {
                N1MC14SUB3_tight2[cbin]->SetBinContent(ebin, N1MCp14SUB3_tight2[cbin]->GetBinContent(ebin));
                N1MC14SUB3_tight2[cbin]->SetBinError(ebin, N1MCp14SUB3_tight2[cbin]->GetBinError(ebin));
            } else {
                N1MC14SUB3_tight2[cbin]->SetBinContent(ebin, N1MCm14SUB3_tight2[cbin]->GetBinContent(ebin));
                N1MC14SUB3_tight2[cbin]->SetBinError(ebin, N1MCm14SUB3_tight2[cbin]->GetBinError(ebin));
            }

        }
    }

    if (!fopen("results/Eta_Distributions","r")) system("mkdir results/Eta_Distributions");

    TH1D * h0 = new TH1D("h0", "", 100, -2.5, 2.5);
    h0->SetStats(0);
    h0->SetXTitle("p_{T} (GeV/c)");
    h0->SetYTitle("v_{1}");

    TCanvas * cN1SUB2_centscan_tight2 = new TCanvas("cN1SUB2_centscan_tight2","cN1SUB2_centscan_tight2",1100,600);
    cN1SUB2_centscan_tight2->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padN1SUB2_centscan_tight2 = (TPad *) cN1SUB2_centscan_tight2->cd(cbin+1);
        padN1SUB2_centscan_tight2->SetGrid();
        TH1D * hN1SUB2_centscan_tight2 = (TH1D *) h0->Clone(Form("hN1SUB2_centscan_tight2_%c",cbin));
        hN1SUB2_centscan_tight2->SetYTitle("v_{1}^{odd}");
        hN1SUB2_centscan_tight2->GetYaxis()->SetRangeUser(-0.012, 0.012);
        hN1SUB2_centscan_tight2->Draw();
        N1SUB2[cbin]->SetMarkerColor(kBlue);
        N1SUB2[cbin]->SetLineColor(kBlue);
        N1SUB2[cbin]->SetMarkerStyle(21);
        N1SUB2[cbin]->SetMarkerSize(1.2);
        N1SUB2[cbin]->Draw("same");
        N1SUB2_tight2[cbin]->SetMarkerColor(kGreen+2);
        N1SUB2_tight2[cbin]->SetLineColor(kGreen+2);
        N1SUB2_tight2[cbin]->SetMarkerStyle(28);
        N1SUB2_tight2[cbin]->SetMarkerSize(1.7);
        //N1SUB2_tight2[cbin]->Draw("same");
        TPaveText * txN1SUB2_centscan_tight2 = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txN1SUB2_centscan_tight2, 18);
        txN1SUB2_centscan_tight2->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txN1SUB2_centscan_tight2->Draw();
    }
    cN1SUB2_centscan_tight2->cd(2);
    TLegend * legN1SUB2_centscan_tight2 = new TLegend(0.05, 0.79, 0.39, 0.97);
    SetLegend(legN1SUB2_centscan_tight2, 16);
    legN1SUB2_centscan_tight2->AddEntry(N1SUB2[0],"nominal","p");
    legN1SUB2_centscan_tight2->AddEntry(N1SUB2_tight2[0],"tight2","p");
    //legN1SUB2_centscan_tight2->Draw();
    cN1SUB2_centscan_tight2->cd(1);
    TPaveText * txN1SUB2_centscan_tight2 = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txN1SUB2_centscan_tight2, 16);
    txN1SUB2_centscan_tight2->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txN1SUB2_centscan_tight2->Draw();
    cN1SUB2_centscan_tight2->Print("results/Eta_Distributions/N1/centscan_N1SUB2_nom_tight2.png","png");

    //
    // TCanvas * cratN1SUB2_centscan_tight2 = new TCanvas("cratN1SUB2_centscan_tight2","cratN1SUB2_centscan_tight2",1100,600);
    // cratN1SUB2_centscan_tight2->Divide(4,2,0,0);
    // for (int cbin = 0; cbin<8; cbin++) {
    //     TPad * padratN1SUB2_centscan_tight2 = (TPad *) cratN1SUB2_centscan_tight2->cd(cbin+1);
    //     padratN1SUB2_centscan_tight2->SetGrid();
    //     TH1D * hratN1SUB2_centscan_tight2 = (TH1D *) h0->Clone(Form("hratN1SUB2_centscan_tight2_%c",cbin));
    //     hratN1SUB2_centscan_tight2->SetYTitle("v_{1}^{odd}{tight2} / v_{1}^{odd}{nominal}");
    //     hratN1SUB2_centscan_tight2->GetYaxis()->SetRangeUser(0.9, 1.1);
    //     hratN1SUB2_centscan_tight2->Draw();
    //     TH1D * ratN1SUB2_centscan_tight2 = (TH1D *) N1SUB2_tight2[cbin]->Clone(Form("ratN1SUB2_%c",cbin));
    //     ratN1SUB2_centscan_tight2->Divide(N1SUB2[cbin]);
    //     for (int k = 1; k<=ratN1SUB2_centscan_tight2->GetNbinsX(); k++) {
    //         double x = N1SUB2_tight2[cbin]->GetBinContent(k);
    //         double y = N1SUB2[cbin]->GetBinContent(k);
    //         double xe = N1SUB2_tight2[cbin]->GetBinError(k);
    //         double ye = N1SUB2[cbin]->GetBinError(k);
    //         ratN1SUB2_centscan_tight2->SetBinError(k, ErrRatCalc(x, y, xe, ye));
    //     }
    //     ratN1SUB2_centscan_tight2->SetMarkerColor(kBlack);
    //     ratN1SUB2_centscan_tight2->SetLineColor(kBlack);
    //     ratN1SUB2_centscan_tight2->SetMarkerStyle(20);
    //     ratN1SUB2_centscan_tight2->SetMarkerSize(1.2);
    //     ratN1SUB2_centscan_tight2->Draw("same");
    //     TPaveText * txratN1SUB2_centscan_tight2 = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
    //     SetTPaveTxt(txratN1SUB2_centscan_tight2, 18);
    //     txratN1SUB2_centscan_tight2->AddText(Form("%d-%d",cminCENT[cbin],cmaxCENT[cbin]));
    //     txratN1SUB2_centscan_tight2->Draw();
    //     TLine * lnratN1SUB2_centscan_tight2 = new TLine(-2.5, 1.0, 2.5, 1.0);
    //     lnratN1SUB2_centscan_tight2->SetLineWidth(2);
    //     lnratN1SUB2_centscan_tight2->Draw();
    // }
    // cratN1SUB2_centscan_tight2->cd(1);
    // TPaveText * txratN1SUB2_centscan_tight2 = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    // SetTPaveTxt(txratN1SUB2_centscan_tight2, 16);
    // txratN1SUB2_centscan_tight2->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    // txratN1SUB2_centscan_tight2->Draw();
    // cratN1SUB2_centscan_tight2->Print("results/Eta_Distributions/N1/centscan_N1SUB2_nom_tight2_ratio.png","png");

}
