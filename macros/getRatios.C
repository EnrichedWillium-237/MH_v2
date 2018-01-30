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
static const int nbinsETA = 14;
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0, -2.4,  0.0};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  0.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24",
                   "-24_0", "0_24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 16;
TString ANAL[] = {"N1MCm22SUB2", "N1MCp22SUB2", "N1MCm22SUB3", "N1MCp22SUB3",
                  "N1SUB2",      "N1ASUB2",     "N1BSUB2",     "N112SUB2",    "N112ASUB2",   "N112BSUB2",
                  "N1SUB3",      "N1ASUB3",     "N1BSUB3"};

# include "src/getRatios.h"

TFile * finPt;
TFile * finEta;
TFile * fout;

void getRatios()
{
    TH1::SetDefaultSumw2();

    finPt = new TFile("hists/MH_combined_Pt.root","read");

    for (int ebin = 0; ebin<nbinsETA; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString mtag = Form("MH_nominal/eta_%s/%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);

            N1MCm22SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCm22SUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            mtag = Form("MH_tight2/eta_%s/%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);

            N1MCm22SUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCm22SUB3_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB3_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB2_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB3_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB3_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB3_tight_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            mtag = Form("MH_wide/eta_%s/%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);

            N1MCm22SUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCm22SUB3_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB3_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB2_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB3_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB3_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB3_wide_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            mtag = Form("MH_narrow/eta_%s/%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);

            N1MCm22SUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCm22SUB3_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB3_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB2_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB3_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB3_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB3_narrow_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
        }
    }

    cout << "Calculating ratios for tight, wide, and narrow cuts... " << endl;
    finEta = new TFile("hists/MH_combined_Eta.root","read");

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TString mtag = Form("MH_nominal/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm22SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        mtag = Form("MH_tight2/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm22SUB3_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB3_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB3_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_tight_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        mtag = Form("MH_wide/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm22SUB3_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB3_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB3_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_wide_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));

        mtag = Form("MH_narrow/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm22SUB3_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB3_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB3_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_narrow_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
    }

    for (int ebin = 0; ebin<nbinsETA; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            rat_N1MCm22SUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N1MCm22SUB2_pT[ebin][cbin]->Clone(Form("rat_N1MCm22SUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCm22SUB2_nom_tight_pT[ebin][cbin]->Divide(N1MCm22SUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCm22SUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCm22SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCm22SUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCm22SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCm22SUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCm22SUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCp22SUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N1MCp22SUB2_pT[ebin][cbin]->Clone(Form("rat_N1MCp22SUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCp22SUB2_nom_tight_pT[ebin][cbin]->Divide(N1MCp22SUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCp22SUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCp22SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCp22SUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCp22SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCp22SUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCp22SUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCm22SUB3_nom_tight_pT[ebin][cbin] = (TH1D *) N1MCm22SUB3_pT[ebin][cbin]->Clone(Form("rat_N1MCm22SUB3_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCm22SUB3_nom_tight_pT[ebin][cbin]->Divide(N1MCm22SUB3_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCm22SUB3_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCm22SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCm22SUB3_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCm22SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCm22SUB3_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCm22SUB3_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCp22SUB3_nom_tight_pT[ebin][cbin] = (TH1D *) N1MCp22SUB3_pT[ebin][cbin]->Clone(Form("rat_N1MCp22SUB3_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCp22SUB3_nom_tight_pT[ebin][cbin]->Divide(N1MCp22SUB3_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCp22SUB3_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCp22SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCp22SUB3_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCp22SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCp22SUB3_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCp22SUB3_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1SUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N1SUB2_pT[ebin][cbin]->Clone(Form("N1SUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1SUB2_nom_tight_pT[ebin][cbin]->Divide(N1SUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1SUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1SUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1SUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1SUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1ASUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N1ASUB2_pT[ebin][cbin]->Clone(Form("N1ASUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1ASUB2_nom_tight_pT[ebin][cbin]->Divide(N1ASUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1ASUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1ASUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1ASUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1ASUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1ASUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1ASUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1BSUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N1BSUB2_pT[ebin][cbin]->Clone(Form("N1BSUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1BSUB2_nom_tight_pT[ebin][cbin]->Divide(N1BSUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1BSUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1BSUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1BSUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1BSUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1BSUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1BSUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112SUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N112SUB2_pT[ebin][cbin]->Clone(Form("N112SUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112SUB2_nom_tight_pT[ebin][cbin]->Divide(N112SUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112SUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112SUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112SUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N112SUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112ASUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N112ASUB2_pT[ebin][cbin]->Clone(Form("N112ASUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112ASUB2_nom_tight_pT[ebin][cbin]->Divide(N112ASUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112ASUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112ASUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112ASUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112ASUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112ASUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N112ASUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112BSUB2_nom_tight_pT[ebin][cbin] = (TH1D *) N112BSUB2_pT[ebin][cbin]->Clone(Form("N112BSUB2_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112BSUB2_nom_tight_pT[ebin][cbin]->Divide(N112BSUB2_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112BSUB2_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112BSUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112BSUB2_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112BSUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112BSUB2_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N112BSUB2_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1SUB3_nom_tight_pT[ebin][cbin] = (TH1D *) N1SUB3_pT[ebin][cbin]->Clone(Form("N1SUB3_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1SUB3_nom_tight_pT[ebin][cbin]->Divide(N1SUB3_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1SUB3_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1SUB3_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1SUB3_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1SUB3_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1ASUB3_nom_tight_pT[ebin][cbin] = (TH1D *) N1ASUB3_pT[ebin][cbin]->Clone(Form("N1ASUB3_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1ASUB3_nom_tight_pT[ebin][cbin]->Divide(N1ASUB3_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1ASUB3_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1ASUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1ASUB3_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1ASUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1ASUB3_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1ASUB3_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1BSUB3_nom_tight_pT[ebin][cbin] = (TH1D *) N1BSUB3_pT[ebin][cbin]->Clone(Form("N1BSUB3_nom_tight_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1BSUB3_nom_tight_pT[ebin][cbin]->Divide(N1BSUB3_tight_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1BSUB3_nom_tight_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1BSUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1BSUB3_tight_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1BSUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1BSUB3_tight_pT[ebin][cbin]->GetBinError(k);
                rat_N1BSUB3_nom_tight_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }

            ///

            rat_N1MCm22SUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N1MCm22SUB2_pT[ebin][cbin]->Clone(Form("rat_N1MCm22SUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCm22SUB2_nom_wide_pT[ebin][cbin]->Divide(N1MCm22SUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCm22SUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCm22SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCm22SUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCm22SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCm22SUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCm22SUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCp22SUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N1MCp22SUB2_pT[ebin][cbin]->Clone(Form("rat_N1MCp22SUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCp22SUB2_nom_wide_pT[ebin][cbin]->Divide(N1MCp22SUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCp22SUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCp22SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCp22SUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCp22SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCp22SUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCp22SUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCm22SUB3_nom_wide_pT[ebin][cbin] = (TH1D *) N1MCm22SUB3_pT[ebin][cbin]->Clone(Form("rat_N1MCm22SUB3_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCm22SUB3_nom_wide_pT[ebin][cbin]->Divide(N1MCm22SUB3_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCm22SUB3_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCm22SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCm22SUB3_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCm22SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCm22SUB3_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCm22SUB3_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCp22SUB3_nom_wide_pT[ebin][cbin] = (TH1D *) N1MCp22SUB3_pT[ebin][cbin]->Clone(Form("rat_N1MCp22SUB3_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCp22SUB3_nom_wide_pT[ebin][cbin]->Divide(N1MCp22SUB3_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCp22SUB3_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCp22SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCp22SUB3_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCp22SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCp22SUB3_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCp22SUB3_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1SUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N1SUB2_pT[ebin][cbin]->Clone(Form("N1SUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1SUB2_nom_wide_pT[ebin][cbin]->Divide(N1SUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1SUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1SUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1SUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1SUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1ASUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N1ASUB2_pT[ebin][cbin]->Clone(Form("N1ASUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1ASUB2_nom_wide_pT[ebin][cbin]->Divide(N1ASUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1ASUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1ASUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1ASUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1ASUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1ASUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1ASUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1BSUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N1BSUB2_pT[ebin][cbin]->Clone(Form("N1BSUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1BSUB2_nom_wide_pT[ebin][cbin]->Divide(N1BSUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1BSUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1BSUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1BSUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1BSUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1BSUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1BSUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112SUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N112SUB2_pT[ebin][cbin]->Clone(Form("N112SUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112SUB2_nom_wide_pT[ebin][cbin]->Divide(N112SUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112SUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112SUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112SUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N112SUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112ASUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N112ASUB2_pT[ebin][cbin]->Clone(Form("N112ASUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112ASUB2_nom_wide_pT[ebin][cbin]->Divide(N112ASUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112ASUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112ASUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112ASUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112ASUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112ASUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N112ASUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112BSUB2_nom_wide_pT[ebin][cbin] = (TH1D *) N112BSUB2_pT[ebin][cbin]->Clone(Form("N112BSUB2_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112BSUB2_nom_wide_pT[ebin][cbin]->Divide(N112BSUB2_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112BSUB2_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112BSUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112BSUB2_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112BSUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112BSUB2_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N112BSUB2_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1SUB3_nom_wide_pT[ebin][cbin] = (TH1D *) N1SUB3_pT[ebin][cbin]->Clone(Form("N1SUB3_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1SUB3_nom_wide_pT[ebin][cbin]->Divide(N1SUB3_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1SUB3_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1SUB3_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1SUB3_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1SUB3_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1ASUB3_nom_wide_pT[ebin][cbin] = (TH1D *) N1ASUB3_pT[ebin][cbin]->Clone(Form("N1ASUB3_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1ASUB3_nom_wide_pT[ebin][cbin]->Divide(N1ASUB3_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1ASUB3_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1ASUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1ASUB3_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1ASUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1ASUB3_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1ASUB3_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1BSUB3_nom_wide_pT[ebin][cbin] = (TH1D *) N1BSUB3_pT[ebin][cbin]->Clone(Form("N1BSUB3_nom_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1BSUB3_nom_wide_pT[ebin][cbin]->Divide(N1BSUB3_wide_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1BSUB3_nom_wide_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1BSUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1BSUB3_wide_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1BSUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1BSUB3_wide_pT[ebin][cbin]->GetBinError(k);
                rat_N1BSUB3_nom_wide_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }

            ///

            rat_N1MCm22SUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N1MCm22SUB2_pT[ebin][cbin]->Clone(Form("rat_N1MCm22SUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCm22SUB2_nom_narrow_pT[ebin][cbin]->Divide(N1MCm22SUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCm22SUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCm22SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCm22SUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCm22SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCm22SUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCm22SUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCp22SUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N1MCp22SUB2_pT[ebin][cbin]->Clone(Form("rat_N1MCp22SUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCp22SUB2_nom_narrow_pT[ebin][cbin]->Divide(N1MCp22SUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCp22SUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCp22SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCp22SUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCp22SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCp22SUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCp22SUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCm22SUB3_nom_narrow_pT[ebin][cbin] = (TH1D *) N1MCm22SUB3_pT[ebin][cbin]->Clone(Form("rat_N1MCm22SUB3_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCm22SUB3_nom_narrow_pT[ebin][cbin]->Divide(N1MCm22SUB3_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCm22SUB3_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCm22SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCm22SUB3_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCm22SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCm22SUB3_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCm22SUB3_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1MCp22SUB3_nom_narrow_pT[ebin][cbin] = (TH1D *) N1MCp22SUB3_pT[ebin][cbin]->Clone(Form("rat_N1MCp22SUB3_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1MCp22SUB3_nom_narrow_pT[ebin][cbin]->Divide(N1MCp22SUB3_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1MCp22SUB3_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1MCp22SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1MCp22SUB3_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1MCp22SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1MCp22SUB3_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1MCp22SUB3_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1SUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N1SUB2_pT[ebin][cbin]->Clone(Form("N1SUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1SUB2_nom_narrow_pT[ebin][cbin]->Divide(N1SUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1SUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1SUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1SUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1SUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1ASUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N1ASUB2_pT[ebin][cbin]->Clone(Form("N1ASUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1ASUB2_nom_narrow_pT[ebin][cbin]->Divide(N1ASUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1ASUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1ASUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1ASUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1ASUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1ASUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1ASUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1BSUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N1BSUB2_pT[ebin][cbin]->Clone(Form("N1BSUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1BSUB2_nom_narrow_pT[ebin][cbin]->Divide(N1BSUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1BSUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1BSUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N1BSUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1BSUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N1BSUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1BSUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112SUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N112SUB2_pT[ebin][cbin]->Clone(Form("N112SUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112SUB2_nom_narrow_pT[ebin][cbin]->Divide(N112SUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112SUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112SUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112SUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112SUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112SUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N112SUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112ASUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N112ASUB2_pT[ebin][cbin]->Clone(Form("N112ASUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112ASUB2_nom_narrow_pT[ebin][cbin]->Divide(N112ASUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112ASUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112ASUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112ASUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112ASUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112ASUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N112ASUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N112BSUB2_nom_narrow_pT[ebin][cbin] = (TH1D *) N112BSUB2_pT[ebin][cbin]->Clone(Form("N112BSUB2_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N112BSUB2_nom_narrow_pT[ebin][cbin]->Divide(N112BSUB2_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N112BSUB2_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N112BSUB2_pT[ebin][cbin]->GetBinContent(k);
                double y = N112BSUB2_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N112BSUB2_pT[ebin][cbin]->GetBinError(k);
                double ye = N112BSUB2_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N112BSUB2_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1SUB3_nom_narrow_pT[ebin][cbin] = (TH1D *) N1SUB3_pT[ebin][cbin]->Clone(Form("N1SUB3_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1SUB3_nom_narrow_pT[ebin][cbin]->Divide(N1SUB3_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1SUB3_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1SUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1SUB3_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1SUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1SUB3_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1SUB3_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1ASUB3_nom_narrow_pT[ebin][cbin] = (TH1D *) N1ASUB3_pT[ebin][cbin]->Clone(Form("N1ASUB3_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1ASUB3_nom_narrow_pT[ebin][cbin]->Divide(N1ASUB3_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1ASUB3_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1ASUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1ASUB3_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1ASUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1ASUB3_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1ASUB3_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
            rat_N1BSUB3_nom_narrow_pT[ebin][cbin] = (TH1D *) N1BSUB3_pT[ebin][cbin]->Clone(Form("N1BSUB3_nom_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            rat_N1BSUB3_nom_narrow_pT[ebin][cbin]->Divide(N1BSUB3_narrow_pT[ebin][cbin]);
            for (int k = 1; k<=rat_N1BSUB3_nom_narrow_pT[ebin][cbin]->GetNbinsX(); k++) {
                double x = N1BSUB3_pT[ebin][cbin]->GetBinContent(k);
                double y = N1BSUB3_narrow_pT[ebin][cbin]->GetBinContent(k);
                double xe = N1BSUB3_pT[ebin][cbin]->GetBinError(k);
                double ye = N1BSUB3_narrow_pT[ebin][cbin]->GetBinError(k);
                rat_N1BSUB3_nom_narrow_pT[ebin][cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
            }
        }
    }

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        rat_N1MCm22SUB2_nom_tight_eta[cbin] = (TH1D *) N1MCm22SUB2_eta[cbin]->Clone(Form("rat_N1MCm22SUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCm22SUB2_nom_tight_eta[cbin]->Divide(N1MCm22SUB2_tight_eta[cbin]);
        for (int k = 0; k<N1MCm22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCm22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MCm22SUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N1MCm22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MCm22SUB2_tight_eta[cbin]->GetBinError(k);
            rat_N1MCm22SUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCp22SUB2_nom_tight_eta[cbin] = (TH1D *) N1MCp22SUB2_eta[cbin]->Clone(Form("rat_N1MCp22SUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCp22SUB2_nom_tight_eta[cbin]->Divide(N1MCp22SUB2_tight_eta[cbin]);
        for (int k = 0; k<N1MCp22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCp22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MCp22SUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N1MCp22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MCp22SUB2_tight_eta[cbin]->GetBinError(k);
            rat_N1MCp22SUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MC22SUB2_nom_tight_eta[cbin] = (TH1D *) N1MC22SUB2_eta[cbin]->Clone(Form("rat_N1MC22SUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MC22SUB2_nom_tight_eta[cbin]->Divide(N1MC22SUB2_tight_eta[cbin]);
        for (int k = 0; k<N1MC22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MC22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MC22SUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N1MC22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MC22SUB2_tight_eta[cbin]->GetBinError(k);
            rat_N1MC22SUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCm22SUB3_nom_tight_eta[cbin] = (TH1D *) N1MCm22SUB3_eta[cbin]->Clone(Form("rat_N1MCm22SUB3_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCm22SUB3_nom_tight_eta[cbin]->Divide(N1MCm22SUB3_tight_eta[cbin]);
        for (int k = 0; k<N1MCm22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCm22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MCm22SUB3_tight_eta[cbin]->GetBinContent(k);
            double xe = N1MCm22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MCm22SUB3_tight_eta[cbin]->GetBinError(k);
            rat_N1MCm22SUB3_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCp22SUB3_nom_tight_eta[cbin] = (TH1D *) N1MCp22SUB3_eta[cbin]->Clone(Form("rat_N1MCp22SUB3_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCp22SUB3_nom_tight_eta[cbin]->Divide(N1MCp22SUB3_tight_eta[cbin]);
        for (int k = 0; k<N1MCp22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCp22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MCp22SUB3_tight_eta[cbin]->GetBinContent(k);
            double xe = N1MCp22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MCp22SUB3_tight_eta[cbin]->GetBinError(k);
            rat_N1MCp22SUB3_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MC22SUB3_nom_tight_eta[cbin] = (TH1D *) N1MC22SUB3_eta[cbin]->Clone(Form("rat_N1MC22SUB3_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MC22SUB3_nom_tight_eta[cbin]->Divide(N1MC22SUB3_tight_eta[cbin]);
        for (int k = 0; k<N1MC22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MC22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MC22SUB3_tight_eta[cbin]->GetBinContent(k);
            double xe = N1MC22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MC22SUB3_tight_eta[cbin]->GetBinError(k);
            rat_N1MC22SUB3_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1SUB2_nom_tight_eta[cbin] = (TH1D *) N1SUB2_eta[cbin]->Clone(Form("rat_N1SUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1SUB2_nom_tight_eta[cbin]->Divide(N1SUB2_tight_eta[cbin]);
        for (int k = 0; k<N1SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1SUB2_eta[cbin]->GetBinContent(k);
            double y = N1SUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N1SUB2_eta[cbin]->GetBinError(k);
            double ye = N1SUB2_tight_eta[cbin]->GetBinError(k);
            rat_N1SUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1ASUB2_nom_tight_eta[cbin] = (TH1D *) N1ASUB2_eta[cbin]->Clone(Form("rat_N1ASUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1ASUB2_nom_tight_eta[cbin]->Divide(N1ASUB2_tight_eta[cbin]);
        for (int k = 0; k<N1ASUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1ASUB2_eta[cbin]->GetBinContent(k);
            double y = N1ASUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N1ASUB2_eta[cbin]->GetBinError(k);
            double ye = N1ASUB2_tight_eta[cbin]->GetBinError(k);
            rat_N1ASUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1BSUB2_nom_tight_eta[cbin] = (TH1D *) N1BSUB2_eta[cbin]->Clone(Form("rat_N1BSUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1BSUB2_nom_tight_eta[cbin]->Divide(N1BSUB2_tight_eta[cbin]);
        for (int k = 0; k<N1BSUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1BSUB2_eta[cbin]->GetBinContent(k);
            double y = N1BSUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N1BSUB2_eta[cbin]->GetBinError(k);
            double ye = N1BSUB2_tight_eta[cbin]->GetBinError(k);
            rat_N1BSUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112SUB2_nom_tight_eta[cbin] = (TH1D *) N112SUB2_eta[cbin]->Clone(Form("rat_N112SUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112SUB2_nom_tight_eta[cbin]->Divide(N112SUB2_tight_eta[cbin]);
        for (int k = 0; k<N112SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112SUB2_eta[cbin]->GetBinContent(k);
            double y = N112SUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N112SUB2_eta[cbin]->GetBinError(k);
            double ye = N112SUB2_tight_eta[cbin]->GetBinError(k);
            rat_N112SUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112ASUB2_nom_tight_eta[cbin] = (TH1D *) N112ASUB2_eta[cbin]->Clone(Form("rat_N112ASUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112ASUB2_nom_tight_eta[cbin]->Divide(N112ASUB2_tight_eta[cbin]);
        for (int k = 0; k<N112ASUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112ASUB2_eta[cbin]->GetBinContent(k);
            double y = N112ASUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N112ASUB2_eta[cbin]->GetBinError(k);
            double ye = N112ASUB2_tight_eta[cbin]->GetBinError(k);
            rat_N112ASUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112BSUB2_nom_tight_eta[cbin] = (TH1D *) N112BSUB2_eta[cbin]->Clone(Form("rat_N112BSUB2_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112BSUB2_nom_tight_eta[cbin]->Divide(N112BSUB2_tight_eta[cbin]);
        for (int k = 0; k<N112BSUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112BSUB2_eta[cbin]->GetBinContent(k);
            double y = N112BSUB2_tight_eta[cbin]->GetBinContent(k);
            double xe = N112BSUB2_eta[cbin]->GetBinError(k);
            double ye = N112BSUB2_tight_eta[cbin]->GetBinError(k);
            rat_N112BSUB2_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1SUB3_nom_tight_eta[cbin] = (TH1D *) N1SUB3_eta[cbin]->Clone(Form("rat_N1SUB3_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1SUB3_nom_tight_eta[cbin]->Divide(N1SUB3_tight_eta[cbin]);
        for (int k = 0; k<N1SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1SUB3_eta[cbin]->GetBinContent(k);
            double y = N1SUB3_tight_eta[cbin]->GetBinContent(k);
            double xe = N1SUB3_eta[cbin]->GetBinError(k);
            double ye = N1SUB3_tight_eta[cbin]->GetBinError(k);
            rat_N1SUB3_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1ASUB3_nom_tight_eta[cbin] = (TH1D *) N1ASUB3_eta[cbin]->Clone(Form("rat_N1ASUB3_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1ASUB3_nom_tight_eta[cbin]->Divide(N1ASUB3_tight_eta[cbin]);
        for (int k = 0; k<N1ASUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1ASUB3_eta[cbin]->GetBinContent(k);
            double y = N1ASUB3_tight_eta[cbin]->GetBinContent(k);
            double xe = N1ASUB3_eta[cbin]->GetBinError(k);
            double ye = N1ASUB3_tight_eta[cbin]->GetBinError(k);
            rat_N1ASUB3_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1BSUB3_nom_tight_eta[cbin] = (TH1D *) N1BSUB3_eta[cbin]->Clone(Form("rat_N1BSUB3_nom_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1BSUB3_nom_tight_eta[cbin]->Divide(N1BSUB3_tight_eta[cbin]);
        for (int k = 0; k<N1BSUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1BSUB3_eta[cbin]->GetBinContent(k);
            double y = N1BSUB3_tight_eta[cbin]->GetBinContent(k);
            double xe = N1BSUB3_eta[cbin]->GetBinError(k);
            double ye = N1BSUB3_tight_eta[cbin]->GetBinError(k);
            rat_N1BSUB3_nom_tight_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }

        ///

        rat_N1MCm22SUB2_nom_wide_eta[cbin] = (TH1D *) N1MCm22SUB2_eta[cbin]->Clone(Form("rat_N1MCm22SUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCm22SUB2_nom_wide_eta[cbin]->Divide(N1MCm22SUB2_wide_eta[cbin]);
        for (int k = 0; k<N1MCm22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCm22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MCm22SUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N1MCm22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MCm22SUB2_wide_eta[cbin]->GetBinError(k);
            rat_N1MCm22SUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCp22SUB2_nom_wide_eta[cbin] = (TH1D *) N1MCp22SUB2_eta[cbin]->Clone(Form("rat_N1MCp22SUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCp22SUB2_nom_wide_eta[cbin]->Divide(N1MCp22SUB2_wide_eta[cbin]);
        for (int k = 0; k<N1MCp22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCp22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MCp22SUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N1MCp22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MCp22SUB2_wide_eta[cbin]->GetBinError(k);
            rat_N1MCp22SUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MC22SUB2_nom_wide_eta[cbin] = (TH1D *) N1MC22SUB2_eta[cbin]->Clone(Form("rat_N1MC22SUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MC22SUB2_nom_wide_eta[cbin]->Divide(N1MC22SUB2_wide_eta[cbin]);
        for (int k = 0; k<N1MC22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MC22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MC22SUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N1MC22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MC22SUB2_wide_eta[cbin]->GetBinError(k);
            rat_N1MC22SUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCm22SUB3_nom_wide_eta[cbin] = (TH1D *) N1MCm22SUB3_eta[cbin]->Clone(Form("rat_N1MCm22SUB3_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCm22SUB3_nom_wide_eta[cbin]->Divide(N1MCm22SUB3_wide_eta[cbin]);
        for (int k = 0; k<N1MCm22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCm22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MCm22SUB3_wide_eta[cbin]->GetBinContent(k);
            double xe = N1MCm22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MCm22SUB3_wide_eta[cbin]->GetBinError(k);
            rat_N1MCm22SUB3_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCp22SUB3_nom_wide_eta[cbin] = (TH1D *) N1MCp22SUB3_eta[cbin]->Clone(Form("rat_N1MCp22SUB3_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCp22SUB3_nom_wide_eta[cbin]->Divide(N1MCp22SUB3_wide_eta[cbin]);
        for (int k = 0; k<N1MCp22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCp22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MCp22SUB3_wide_eta[cbin]->GetBinContent(k);
            double xe = N1MCp22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MCp22SUB3_wide_eta[cbin]->GetBinError(k);
            rat_N1MCp22SUB3_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MC22SUB3_nom_wide_eta[cbin] = (TH1D *) N1MC22SUB3_eta[cbin]->Clone(Form("rat_N1MC22SUB3_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MC22SUB3_nom_wide_eta[cbin]->Divide(N1MC22SUB3_wide_eta[cbin]);
        for (int k = 0; k<N1MC22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MC22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MC22SUB3_wide_eta[cbin]->GetBinContent(k);
            double xe = N1MC22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MC22SUB3_wide_eta[cbin]->GetBinError(k);
            rat_N1MC22SUB3_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1SUB2_nom_wide_eta[cbin] = (TH1D *) N1SUB2_eta[cbin]->Clone(Form("rat_N1SUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1SUB2_nom_wide_eta[cbin]->Divide(N1SUB2_wide_eta[cbin]);
        for (int k = 0; k<N1SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1SUB2_eta[cbin]->GetBinContent(k);
            double y = N1SUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N1SUB2_eta[cbin]->GetBinError(k);
            double ye = N1SUB2_wide_eta[cbin]->GetBinError(k);
            rat_N1SUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1ASUB2_nom_wide_eta[cbin] = (TH1D *) N1ASUB2_eta[cbin]->Clone(Form("rat_N1ASUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1ASUB2_nom_wide_eta[cbin]->Divide(N1ASUB2_wide_eta[cbin]);
        for (int k = 0; k<N1ASUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1ASUB2_eta[cbin]->GetBinContent(k);
            double y = N1ASUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N1ASUB2_eta[cbin]->GetBinError(k);
            double ye = N1ASUB2_wide_eta[cbin]->GetBinError(k);
            rat_N1ASUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1BSUB2_nom_wide_eta[cbin] = (TH1D *) N1BSUB2_eta[cbin]->Clone(Form("rat_N1BSUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1BSUB2_nom_wide_eta[cbin]->Divide(N1BSUB2_wide_eta[cbin]);
        for (int k = 0; k<N1BSUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1BSUB2_eta[cbin]->GetBinContent(k);
            double y = N1BSUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N1BSUB2_eta[cbin]->GetBinError(k);
            double ye = N1BSUB2_wide_eta[cbin]->GetBinError(k);
            rat_N1BSUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112SUB2_nom_wide_eta[cbin] = (TH1D *) N112SUB2_eta[cbin]->Clone(Form("rat_N112SUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112SUB2_nom_wide_eta[cbin]->Divide(N112SUB2_wide_eta[cbin]);
        for (int k = 0; k<N112SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112SUB2_eta[cbin]->GetBinContent(k);
            double y = N112SUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N112SUB2_eta[cbin]->GetBinError(k);
            double ye = N112SUB2_wide_eta[cbin]->GetBinError(k);
            rat_N112SUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112ASUB2_nom_wide_eta[cbin] = (TH1D *) N112ASUB2_eta[cbin]->Clone(Form("rat_N112ASUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112ASUB2_nom_wide_eta[cbin]->Divide(N112ASUB2_wide_eta[cbin]);
        for (int k = 0; k<N112ASUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112ASUB2_eta[cbin]->GetBinContent(k);
            double y = N112ASUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N112ASUB2_eta[cbin]->GetBinError(k);
            double ye = N112ASUB2_wide_eta[cbin]->GetBinError(k);
            rat_N112ASUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112BSUB2_nom_wide_eta[cbin] = (TH1D *) N112BSUB2_eta[cbin]->Clone(Form("rat_N112BSUB2_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112BSUB2_nom_wide_eta[cbin]->Divide(N112BSUB2_wide_eta[cbin]);
        for (int k = 0; k<N112BSUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112BSUB2_eta[cbin]->GetBinContent(k);
            double y = N112BSUB2_wide_eta[cbin]->GetBinContent(k);
            double xe = N112BSUB2_eta[cbin]->GetBinError(k);
            double ye = N112BSUB2_wide_eta[cbin]->GetBinError(k);
            rat_N112BSUB2_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1SUB3_nom_wide_eta[cbin] = (TH1D *) N1SUB3_eta[cbin]->Clone(Form("rat_N1SUB3_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1SUB3_nom_wide_eta[cbin]->Divide(N1SUB3_wide_eta[cbin]);
        for (int k = 0; k<N1SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1SUB3_eta[cbin]->GetBinContent(k);
            double y = N1SUB3_wide_eta[cbin]->GetBinContent(k);
            double xe = N1SUB3_eta[cbin]->GetBinError(k);
            double ye = N1SUB3_wide_eta[cbin]->GetBinError(k);
            rat_N1SUB3_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1ASUB3_nom_wide_eta[cbin] = (TH1D *) N1ASUB3_eta[cbin]->Clone(Form("rat_N1ASUB3_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1ASUB3_nom_wide_eta[cbin]->Divide(N1ASUB3_wide_eta[cbin]);
        for (int k = 0; k<N1ASUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1ASUB3_eta[cbin]->GetBinContent(k);
            double y = N1ASUB3_wide_eta[cbin]->GetBinContent(k);
            double xe = N1ASUB3_eta[cbin]->GetBinError(k);
            double ye = N1ASUB3_wide_eta[cbin]->GetBinError(k);
            rat_N1ASUB3_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1BSUB3_nom_wide_eta[cbin] = (TH1D *) N1BSUB3_eta[cbin]->Clone(Form("rat_N1BSUB3_nom_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1BSUB3_nom_wide_eta[cbin]->Divide(N1BSUB3_wide_eta[cbin]);
        for (int k = 0; k<N1BSUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1BSUB3_eta[cbin]->GetBinContent(k);
            double y = N1BSUB3_wide_eta[cbin]->GetBinContent(k);
            double xe = N1BSUB3_eta[cbin]->GetBinError(k);
            double ye = N1BSUB3_wide_eta[cbin]->GetBinError(k);
            rat_N1BSUB3_nom_wide_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }

        ///

        rat_N1MCm22SUB2_nom_narrow_eta[cbin] = (TH1D *) N1MCm22SUB2_eta[cbin]->Clone(Form("rat_N1MCm22SUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCm22SUB2_nom_narrow_eta[cbin]->Divide(N1MCm22SUB2_narrow_eta[cbin]);
        for (int k = 0; k<N1MCm22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCm22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MCm22SUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1MCm22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MCm22SUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N1MCm22SUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCp22SUB2_nom_narrow_eta[cbin] = (TH1D *) N1MCp22SUB2_eta[cbin]->Clone(Form("rat_N1MCp22SUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCp22SUB2_nom_narrow_eta[cbin]->Divide(N1MCp22SUB2_narrow_eta[cbin]);
        for (int k = 0; k<N1MCp22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCp22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MCp22SUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1MCp22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MCp22SUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N1MCp22SUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MC22SUB2_nom_narrow_eta[cbin] = (TH1D *) N1MC22SUB2_eta[cbin]->Clone(Form("rat_N1MC22SUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MC22SUB2_nom_narrow_eta[cbin]->Divide(N1MC22SUB2_narrow_eta[cbin]);
        for (int k = 0; k<N1MC22SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MC22SUB2_eta[cbin]->GetBinContent(k);
            double y = N1MC22SUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1MC22SUB2_eta[cbin]->GetBinError(k);
            double ye = N1MC22SUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N1MC22SUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCm22SUB3_nom_narrow_eta[cbin] = (TH1D *) N1MCm22SUB3_eta[cbin]->Clone(Form("rat_N1MCm22SUB3_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCm22SUB3_nom_narrow_eta[cbin]->Divide(N1MCm22SUB3_narrow_eta[cbin]);
        for (int k = 0; k<N1MCm22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCm22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MCm22SUB3_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1MCm22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MCm22SUB3_narrow_eta[cbin]->GetBinError(k);
            rat_N1MCm22SUB3_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MCp22SUB3_nom_narrow_eta[cbin] = (TH1D *) N1MCp22SUB3_eta[cbin]->Clone(Form("rat_N1MCp22SUB3_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MCp22SUB3_nom_narrow_eta[cbin]->Divide(N1MCp22SUB3_narrow_eta[cbin]);
        for (int k = 0; k<N1MCp22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MCp22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MCp22SUB3_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1MCp22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MCp22SUB3_narrow_eta[cbin]->GetBinError(k);
            rat_N1MCp22SUB3_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1MC22SUB3_nom_narrow_eta[cbin] = (TH1D *) N1MC22SUB3_eta[cbin]->Clone(Form("rat_N1MC22SUB3_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1MC22SUB3_nom_narrow_eta[cbin]->Divide(N1MC22SUB3_narrow_eta[cbin]);
        for (int k = 0; k<N1MC22SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1MC22SUB3_eta[cbin]->GetBinContent(k);
            double y = N1MC22SUB3_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1MC22SUB3_eta[cbin]->GetBinError(k);
            double ye = N1MC22SUB3_narrow_eta[cbin]->GetBinError(k);
            rat_N1MC22SUB3_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1SUB2_nom_narrow_eta[cbin] = (TH1D *) N1SUB2_eta[cbin]->Clone(Form("rat_N1SUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1SUB2_nom_narrow_eta[cbin]->Divide(N1SUB2_narrow_eta[cbin]);
        for (int k = 0; k<N1SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1SUB2_eta[cbin]->GetBinContent(k);
            double y = N1SUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1SUB2_eta[cbin]->GetBinError(k);
            double ye = N1SUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N1SUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1ASUB2_nom_narrow_eta[cbin] = (TH1D *) N1ASUB2_eta[cbin]->Clone(Form("rat_N1ASUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1ASUB2_nom_narrow_eta[cbin]->Divide(N1ASUB2_narrow_eta[cbin]);
        for (int k = 0; k<N1ASUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1ASUB2_eta[cbin]->GetBinContent(k);
            double y = N1ASUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1ASUB2_eta[cbin]->GetBinError(k);
            double ye = N1ASUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N1ASUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1BSUB2_nom_narrow_eta[cbin] = (TH1D *) N1BSUB2_eta[cbin]->Clone(Form("rat_N1BSUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1BSUB2_nom_narrow_eta[cbin]->Divide(N1BSUB2_narrow_eta[cbin]);
        for (int k = 0; k<N1BSUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N1BSUB2_eta[cbin]->GetBinContent(k);
            double y = N1BSUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1BSUB2_eta[cbin]->GetBinError(k);
            double ye = N1BSUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N1BSUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112SUB2_nom_narrow_eta[cbin] = (TH1D *) N112SUB2_eta[cbin]->Clone(Form("rat_N112SUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112SUB2_nom_narrow_eta[cbin]->Divide(N112SUB2_narrow_eta[cbin]);
        for (int k = 0; k<N112SUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112SUB2_eta[cbin]->GetBinContent(k);
            double y = N112SUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N112SUB2_eta[cbin]->GetBinError(k);
            double ye = N112SUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N112SUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112ASUB2_nom_narrow_eta[cbin] = (TH1D *) N112ASUB2_eta[cbin]->Clone(Form("rat_N112ASUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112ASUB2_nom_narrow_eta[cbin]->Divide(N112ASUB2_narrow_eta[cbin]);
        for (int k = 0; k<N112ASUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112ASUB2_eta[cbin]->GetBinContent(k);
            double y = N112ASUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N112ASUB2_eta[cbin]->GetBinError(k);
            double ye = N112ASUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N112ASUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N112BSUB2_nom_narrow_eta[cbin] = (TH1D *) N112BSUB2_eta[cbin]->Clone(Form("rat_N112BSUB2_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N112BSUB2_nom_narrow_eta[cbin]->Divide(N112BSUB2_narrow_eta[cbin]);
        for (int k = 0; k<N112BSUB2_eta[cbin]->GetNbinsX(); k++) {
            double x = N112BSUB2_eta[cbin]->GetBinContent(k);
            double y = N112BSUB2_narrow_eta[cbin]->GetBinContent(k);
            double xe = N112BSUB2_eta[cbin]->GetBinError(k);
            double ye = N112BSUB2_narrow_eta[cbin]->GetBinError(k);
            rat_N112BSUB2_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1SUB3_nom_narrow_eta[cbin] = (TH1D *) N1SUB3_eta[cbin]->Clone(Form("rat_N1SUB3_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1SUB3_nom_narrow_eta[cbin]->Divide(N1SUB3_narrow_eta[cbin]);
        for (int k = 0; k<N1SUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1SUB3_eta[cbin]->GetBinContent(k);
            double y = N1SUB3_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1SUB3_eta[cbin]->GetBinError(k);
            double ye = N1SUB3_narrow_eta[cbin]->GetBinError(k);
            rat_N1SUB3_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1ASUB3_nom_narrow_eta[cbin] = (TH1D *) N1ASUB3_eta[cbin]->Clone(Form("rat_N1ASUB3_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1ASUB3_nom_narrow_eta[cbin]->Divide(N1ASUB3_narrow_eta[cbin]);
        for (int k = 0; k<N1ASUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1ASUB3_eta[cbin]->GetBinContent(k);
            double y = N1ASUB3_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1ASUB3_eta[cbin]->GetBinError(k);
            double ye = N1ASUB3_narrow_eta[cbin]->GetBinError(k);
            rat_N1ASUB3_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
        rat_N1BSUB3_nom_narrow_eta[cbin] = (TH1D *) N1BSUB3_eta[cbin]->Clone(Form("rat_N1BSUB3_nom_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        rat_N1BSUB3_nom_narrow_eta[cbin]->Divide(N1BSUB3_narrow_eta[cbin]);
        for (int k = 0; k<N1BSUB3_eta[cbin]->GetNbinsX(); k++) {
            double x = N1BSUB3_eta[cbin]->GetBinContent(k);
            double y = N1BSUB3_narrow_eta[cbin]->GetBinContent(k);
            double xe = N1BSUB3_eta[cbin]->GetBinError(k);
            double ye = N1BSUB3_narrow_eta[cbin]->GetBinError(k);
            rat_N1BSUB3_nom_narrow_eta[cbin]->SetBinError(k, ErrRatCalc( x, y, xe, ye ));
        }
    }


    fout = new TFile("hists/MH_ratios.root","recreate");
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TDirectory * tdcent = (TDirectory *) fout->mkdir(Form("%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        TDirectory * tdPt = (TDirectory *) tdcent->mkdir("vn_pT");
        for (int ebin = 0; ebin<nbinsETA; ebin++) {
            TDirectory * tdPt_ebin = (TDirectory *) tdPt->mkdir(Form("eta_%s",etags[ebin].Data()));
            tdPt_ebin->cd();

            rat_N1MCm22SUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N1MCp22SUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N1MCm22SUB3_nom_tight_pT[ebin][cbin]->Write();
            rat_N1MCp22SUB3_nom_tight_pT[ebin][cbin]->Write();
            rat_N1SUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N1ASUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N1BSUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N112SUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N112ASUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N112BSUB2_nom_tight_pT[ebin][cbin]->Write();
            rat_N1SUB3_nom_tight_pT[ebin][cbin]->Write();
            rat_N1ASUB3_nom_tight_pT[ebin][cbin]->Write();
            rat_N1BSUB3_nom_tight_pT[ebin][cbin]->Write();

            rat_N1MCm22SUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N1MCp22SUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N1MCm22SUB3_nom_wide_pT[ebin][cbin]->Write();
            rat_N1MCp22SUB3_nom_wide_pT[ebin][cbin]->Write();
            rat_N1SUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N1ASUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N1BSUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N112SUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N112ASUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N112BSUB2_nom_wide_pT[ebin][cbin]->Write();
            rat_N1SUB3_nom_wide_pT[ebin][cbin]->Write();
            rat_N1ASUB3_nom_wide_pT[ebin][cbin]->Write();
            rat_N1BSUB3_nom_wide_pT[ebin][cbin]->Write();

            rat_N1MCm22SUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1MCp22SUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1MCm22SUB3_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1MCp22SUB3_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1SUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1ASUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1BSUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N112SUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N112ASUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N112BSUB2_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1SUB3_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1ASUB3_nom_narrow_pT[ebin][cbin]->Write();
            rat_N1BSUB3_nom_narrow_pT[ebin][cbin]->Write();
        }
        TDirectory * tdEta = (TDirectory *) tdcent->mkdir("vn_eta");
        tdEta->cd();

        rat_N1MCm22SUB2_nom_tight_eta[cbin]->Write();
        rat_N1MCp22SUB2_nom_tight_eta[cbin]->Write();
        rat_N1MC22SUB2_nom_tight_eta[cbin]->Write();
        rat_N1MCm22SUB3_nom_tight_eta[cbin]->Write();
        rat_N1MCp22SUB3_nom_tight_eta[cbin]->Write();
        rat_N1MC22SUB3_nom_tight_eta[cbin]->Write();
        rat_N1SUB2_nom_tight_eta[cbin]->Write();
        rat_N1ASUB2_nom_tight_eta[cbin]->Write();
        rat_N1BSUB2_nom_tight_eta[cbin]->Write();
        rat_N112SUB2_nom_tight_eta[cbin]->Write();
        rat_N112ASUB2_nom_tight_eta[cbin]->Write();
        rat_N112BSUB2_nom_tight_eta[cbin]->Write();
        rat_N1SUB3_nom_tight_eta[cbin]->Write();
        rat_N1ASUB3_nom_tight_eta[cbin]->Write();
        rat_N1BSUB3_nom_tight_eta[cbin]->Write();

        rat_N1MCm22SUB2_nom_wide_eta[cbin]->Write();
        rat_N1MCp22SUB2_nom_wide_eta[cbin]->Write();
        rat_N1MC22SUB2_nom_wide_eta[cbin]->Write();
        rat_N1MCm22SUB3_nom_wide_eta[cbin]->Write();
        rat_N1MCp22SUB3_nom_wide_eta[cbin]->Write();
        rat_N1MC22SUB3_nom_wide_eta[cbin]->Write();
        rat_N1SUB2_nom_wide_eta[cbin]->Write();
        rat_N1ASUB2_nom_wide_eta[cbin]->Write();
        rat_N1BSUB2_nom_wide_eta[cbin]->Write();
        rat_N112SUB2_nom_wide_eta[cbin]->Write();
        rat_N112ASUB2_nom_wide_eta[cbin]->Write();
        rat_N112BSUB2_nom_wide_eta[cbin]->Write();
        rat_N1SUB3_nom_wide_eta[cbin]->Write();
        rat_N1ASUB3_nom_wide_eta[cbin]->Write();
        rat_N1BSUB3_nom_wide_eta[cbin]->Write();

        rat_N1MCm22SUB2_nom_narrow_eta[cbin]->Write();
        rat_N1MCp22SUB2_nom_narrow_eta[cbin]->Write();
        rat_N1MC22SUB2_nom_narrow_eta[cbin]->Write();
        rat_N1MCm22SUB3_nom_narrow_eta[cbin]->Write();
        rat_N1MCp22SUB3_nom_narrow_eta[cbin]->Write();
        rat_N1MC22SUB3_nom_narrow_eta[cbin]->Write();
        rat_N1SUB2_nom_narrow_eta[cbin]->Write();
        rat_N1ASUB2_nom_narrow_eta[cbin]->Write();
        rat_N1BSUB2_nom_narrow_eta[cbin]->Write();
        rat_N112SUB2_nom_narrow_eta[cbin]->Write();
        rat_N112ASUB2_nom_narrow_eta[cbin]->Write();
        rat_N112BSUB2_nom_narrow_eta[cbin]->Write();
        rat_N1SUB3_nom_narrow_eta[cbin]->Write();
        rat_N1ASUB3_nom_narrow_eta[cbin]->Write();
        rat_N1BSUB3_nom_narrow_eta[cbin]->Write();
    }

    fout->Close();

    cout << "vn cut ratios written out to hists/MH_ratios.root \n" << endl;

}
