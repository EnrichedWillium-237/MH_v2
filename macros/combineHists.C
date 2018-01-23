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
// static const int netabins = 12;
// static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const int netabins = 14;
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0, -2.4,  0.0};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  0.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24", "-24_0", "0_-24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 12;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2",
                  "N1MCp22SUB2", "N1MCp18SUB2", "N1MCp14SUB2",
                  "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3",
                  "N1MCp22SUB3", "N1MCp18SUB3", "N1MCp14SUB3",
};

TGraphErrors * gN1MCm22SUB2[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB2[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB2[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB2[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB2[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB2[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB3[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB3[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB3[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB3[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB3[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB3[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB2_tight2[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB2_tight2[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB2_tight2[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB2_tight2[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB2_tight2[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB2_tight2[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB3_tight2[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB3_tight2[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB3_tight2[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB3_tight2[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB3_tight2[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB3_tight2[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB2_narrow[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB2_narrow[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB2_narrow[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB2_narrow[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB2_narrow[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB2_narrow[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB3_narrow[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB3_narrow[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB3_narrow[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB3_narrow[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB3_narrow[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB3_narrow[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB2_wide[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB2_wide[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB2_wide[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB2_wide[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB2_wide[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB2_wide[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB3_wide[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB3_wide[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB3_wide[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB3_wide[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB3_wide[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB3_wide[netabins][ncentbins];

TH1D * hN1MCm22SUB2[netabins][ncentbins];
TH1D * hN1MCm18SUB2[netabins][ncentbins];
TH1D * hN1MCm14SUB2[netabins][ncentbins];
TH1D * hN1MCp22SUB2[netabins][ncentbins];
TH1D * hN1MCp18SUB2[netabins][ncentbins];
TH1D * hN1MCp14SUB2[netabins][ncentbins];

TH1D * hN1MCm22SUB3[netabins][ncentbins];
TH1D * hN1MCm18SUB3[netabins][ncentbins];
TH1D * hN1MCm14SUB3[netabins][ncentbins];
TH1D * hN1MCp22SUB3[netabins][ncentbins];
TH1D * hN1MCp18SUB3[netabins][ncentbins];
TH1D * hN1MCp14SUB3[netabins][ncentbins];

TH1D * hN1MCm22SUB2_tight2[netabins][ncentbins];
TH1D * hN1MCm18SUB2_tight2[netabins][ncentbins];
TH1D * hN1MCm14SUB2_tight2[netabins][ncentbins];
TH1D * hN1MCp22SUB2_tight2[netabins][ncentbins];
TH1D * hN1MCp18SUB2_tight2[netabins][ncentbins];
TH1D * hN1MCp14SUB2_tight2[netabins][ncentbins];

TH1D * hN1MCm22SUB3_tight2[netabins][ncentbins];
TH1D * hN1MCm18SUB3_tight2[netabins][ncentbins];
TH1D * hN1MCm14SUB3_tight2[netabins][ncentbins];
TH1D * hN1MCp22SUB3_tight2[netabins][ncentbins];
TH1D * hN1MCp18SUB3_tight2[netabins][ncentbins];
TH1D * hN1MCp14SUB3_tight2[netabins][ncentbins];

TH1D * hN1MCm22SUB2_narrow[netabins][ncentbins];
TH1D * hN1MCm18SUB2_narrow[netabins][ncentbins];
TH1D * hN1MCm14SUB2_narrow[netabins][ncentbins];
TH1D * hN1MCp22SUB2_narrow[netabins][ncentbins];
TH1D * hN1MCp18SUB2_narrow[netabins][ncentbins];
TH1D * hN1MCp14SUB2_narrow[netabins][ncentbins];

TH1D * hN1MCm22SUB3_narrow[netabins][ncentbins];
TH1D * hN1MCm18SUB3_narrow[netabins][ncentbins];
TH1D * hN1MCm14SUB3_narrow[netabins][ncentbins];
TH1D * hN1MCp22SUB3_narrow[netabins][ncentbins];
TH1D * hN1MCp18SUB3_narrow[netabins][ncentbins];
TH1D * hN1MCp14SUB3_narrow[netabins][ncentbins];

TH1D * hN1MCm22SUB2_wide[netabins][ncentbins];
TH1D * hN1MCm18SUB2_wide[netabins][ncentbins];
TH1D * hN1MCm14SUB2_wide[netabins][ncentbins];
TH1D * hN1MCp22SUB2_wide[netabins][ncentbins];
TH1D * hN1MCp18SUB2_wide[netabins][ncentbins];
TH1D * hN1MCp14SUB2_wide[netabins][ncentbins];

TH1D * hN1MCm22SUB3_wide[netabins][ncentbins];
TH1D * hN1MCm18SUB3_wide[netabins][ncentbins];
TH1D * hN1MCm14SUB3_wide[netabins][ncentbins];
TH1D * hN1MCp22SUB3_wide[netabins][ncentbins];
TH1D * hN1MCp18SUB3_wide[netabins][ncentbins];
TH1D * hN1MCp14SUB3_wide[netabins][ncentbins];

TFile * tfin_nominal;
TFile * tfin_tight2;
TFile * tfin_narrow;
TFile * tfin_wide;

TFile * tfoutPt;
TFile * tfoutEta;

void GetPtDists()
{
    tfoutPt = new TFile("results/MH_combined_Pt.root","recreate");

    for (int i = 0; i<NANALS; i++) {
        for (int ebin = 0; ebin<netabins; ebin++) {
            for (int cbin = 0; cbin<ncentbins; cbin++) {
                TString mtag = Form("%s/%s/%d_%d",ANAL[i].Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);
                //cout<<mtag.Data()<<endl;
                if (i == 0) gN1MCm22SUB2[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 1) gN1MCm18SUB2[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 2) gN1MCm14SUB2[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 3) gN1MCp22SUB2[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 4) gN1MCp18SUB2[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 5) gN1MCp14SUB2[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));

                if (i == 6) gN1MCm22SUB3[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 7) gN1MCm18SUB3[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 8) gN1MCm14SUB3[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 9) gN1MCp22SUB3[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 10) gN1MCp18SUB3[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));
                if (i == 11) gN1MCp14SUB3[cbin][cbin] = (TGraphErrors *) tfin_nominal->Get(Form("%s/gA",mtag.data()));


                if (i == 0) gN1MCm22SUB2_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 1) gN1MCm18SUB2_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 2) gN1MCm14SUB2_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 3) gN1MCp22SUB2_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 4) gN1MCp18SUB2_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 5) gN1MCp14SUB2_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));

                if (i == 6) gN1MCm22SUB3_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 7) gN1MCm18SUB3_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 8) gN1MCm14SUB3_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 9) gN1MCp22SUB3_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 10) gN1MCp18SUB3_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));
                if (i == 11) gN1MCp14SUB3_tight2[cbin][cbin] = (TGraphErrors *) tfin_tight2->Get(Form("%s/gA",mtag.data()));


                if (i == 0) gN1MCm22SUB2_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 1) gN1MCm18SUB2_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 2) gN1MCm14SUB2_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 3) gN1MCp22SUB2_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 4) gN1MCp18SUB2_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 5) gN1MCp14SUB2_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));

                if (i == 6) gN1MCm22SUB3_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 7) gN1MCm18SUB3_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 8) gN1MCm14SUB3_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 9) gN1MCp22SUB3_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 10) gN1MCp18SUB3_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));
                if (i == 11) gN1MCp14SUB3_narrow[cbin][cbin] = (TGraphErrors *) tfin_narrow->Get(Form("%s/gA",mtag.data()));


                if (i == 0) gN1MCm22SUB2_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 1) gN1MCm18SUB2_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 2) gN1MCm14SUB2_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 3) gN1MCp22SUB2_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 4) gN1MCp18SUB2_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 5) gN1MCp14SUB2_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));

                if (i == 6) gN1MCm22SUB3_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 7) gN1MCm18SUB3_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 8) gN1MCm14SUB3_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 9) gN1MCp22SUB3_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 10) gN1MCp18SUB3_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
                if (i == 11) gN1MCp14SUB3_wide[cbin][cbin] = (TGraphErrors *) tfin_wide->Get(Form("%s/gA",mtag.data()));
            }
        }
    }

    for (int ebin = 0; ebin<netabins; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            hN1MCm22SUB2[ebin][cbin] = new TH1D(Form("N1MCm22SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            hN1MCm22SUB2[ebin][cbin]->SetStats(0);
            hN1MCm22SUB2[ebin][cbin]->SetXTitle("p_{T} (GeV/c)");
            hN1MCm22SUB2[ebin][cbin]->SetYTitle("v_{1}^{even}");

            hN1MCm18SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm18SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm14SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp22SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp18SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp14SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            hN1MCm18SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm18SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm14SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp22SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp18SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp14SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));


            hN1MCm22SUB2_tight2[ebin][cbin] = new TH1D(Form("N1MCm22SUB2_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            hN1MCm22SUB2_tight2[ebin][cbin]->SetStats(0);
            hN1MCm22SUB2_tight2[ebin][cbin]->SetXTitle("p_{T} (GeV/c)");
            hN1MCm22SUB2_tight2[ebin][cbin]->SetYTitle("v_{1}^{even} (tight2 cuts)");

            hN1MCm18SUB2_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCm18SUB2_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB2_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCm14SUB2_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB2_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCp22SUB2_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB2_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCp18SUB2_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB2_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCp14SUB2_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            hN1MCm18SUB3_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCm18SUB3_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB3_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCm14SUB3_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB3_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCp22SUB3_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB3_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCp18SUB3_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB3_tight2[ebin][cbin] = (TH1D *) hN1MCm22SUB2_tight2[ebin][cbin]->Clone(Form("N1MCp14SUB3_tight2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));


            hN1MCm22SUB2_narrow[ebin][cbin] = new TH1D(Form("N1MCm22SUB2_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            hN1MCm22SUB2_narrow[ebin][cbin]->SetStats(0);
            hN1MCm22SUB2_narrow[ebin][cbin]->SetXTitle("p_{T} (GeV/c)");
            hN1MCm22SUB2_narrow[ebin][cbin]->SetYTitle("v_{1}^{even} (narrow cuts)");

            hN1MCm18SUB2_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCm18SUB2_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB2_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCm14SUB2_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB2_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCp22SUB2_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB2_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCp18SUB2_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB2_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCp14SUB2_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            hN1MCm18SUB3_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCm18SUB3_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB3_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCm14SUB3_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB3_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCp22SUB3_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB3_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCp18SUB3_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB3_narrow[ebin][cbin] = (TH1D *) hN1MCm22SUB2_narrow[ebin][cbin]->Clone(Form("N1MCp14SUB3_narrow_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));


            hN1MCm22SUB2_wide[ebin][cbin] = new TH1D(Form("N1MCm22SUB2_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            hN1MCm22SUB2_wide[ebin][cbin]->SetStats(0);
            hN1MCm22SUB2_wide[ebin][cbin]->SetXTitle("p_{T} (GeV/c)");
            hN1MCm22SUB2_wide[ebin][cbin]->SetYTitle("v_{1}^{even} (wide cuts)");

            hN1MCm18SUB2_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCm18SUB2_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB2_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCm14SUB2_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB2_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCp22SUB2_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB2_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCp18SUB2_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB2_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCp14SUB2_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            hN1MCm18SUB3_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCm18SUB3_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB3_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCm14SUB3_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB3_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCp22SUB3_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB3_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCp18SUB3_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB3_wide[ebin][cbin] = (TH1D *) hN1MCm22SUB2_wide[ebin][cbin]->Clone(Form("N1MCp14SUB3_wide_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));


            GraphToHist( gN1MCm22SUB2[ebin][cbin], hN1MCm22SUB2[ebin][cbin] );
            GraphToHist( gN1MCm18SUB2[ebin][cbin], hN1MCp18SUB2[ebin][cbin] );
            GraphToHist( gN1MCm14SUB2[ebin][cbin], hN1MCp14SUB2[ebin][cbin] );
            GraphToHist( gN1MCp22SUB2[ebin][cbin], hN1MCp22SUB2[ebin][cbin] );
            GraphToHist( gN1MCp18SUB2[ebin][cbin], hN1MCp18SUB2[ebin][cbin] );
            GraphToHist( gN1MCp14SUB2[ebin][cbin], hN1MCp14SUB2[ebin][cbin] );

            GraphToHist( gN1MCm22SUB3[ebin][cbin], hN1MCm22SUB3[ebin][cbin] );
            GraphToHist( gN1MCm18SUB3[ebin][cbin], hN1MCp18SUB3[ebin][cbin] );
            GraphToHist( gN1MCm14SUB3[ebin][cbin], hN1MCp14SUB3[ebin][cbin] );
            GraphToHist( gN1MCp22SUB3[ebin][cbin], hN1MCp22SUB3[ebin][cbin] );
            GraphToHist( gN1MCp18SUB3[ebin][cbin], hN1MCp18SUB3[ebin][cbin] );
            GraphToHist( gN1MCp14SUB3[ebin][cbin], hN1MCp14SUB3[ebin][cbin] );


            GraphToHist( gN1MCm22SUB2_tight2[ebin][cbin], hN1MCm22SUB2_tight2[ebin][cbin] );
            GraphToHist( gN1MCm18SUB2_tight2[ebin][cbin], hN1MCp18SUB2_tight2[ebin][cbin] );
            GraphToHist( gN1MCm14SUB2_tight2[ebin][cbin], hN1MCp14SUB2_tight2[ebin][cbin] );
            GraphToHist( gN1MCp22SUB2_tight2[ebin][cbin], hN1MCp22SUB2_tight2[ebin][cbin] );
            GraphToHist( gN1MCp18SUB2_tight2[ebin][cbin], hN1MCp18SUB2_tight2[ebin][cbin] );
            GraphToHist( gN1MCp14SUB2_tight2[ebin][cbin], hN1MCp14SUB2_tight2[ebin][cbin] );

            GraphToHist( gN1MCm22SUB3_tight2[ebin][cbin], hN1MCm22SUB3_tight2[ebin][cbin] );
            GraphToHist( gN1MCm18SUB3_tight2[ebin][cbin], hN1MCp18SUB3_tight2[ebin][cbin] );
            GraphToHist( gN1MCm14SUB3_tight2[ebin][cbin], hN1MCp14SUB3_tight2[ebin][cbin] );
            GraphToHist( gN1MCp22SUB3_tight2[ebin][cbin], hN1MCp22SUB3_tight2[ebin][cbin] );
            GraphToHist( gN1MCp18SUB3_tight2[ebin][cbin], hN1MCp18SUB3_tight2[ebin][cbin] );
            GraphToHist( gN1MCp14SUB3_tight2[ebin][cbin], hN1MCp14SUB3_tight2[ebin][cbin] );


            GraphToHist( gN1MCm22SUB2_narrow[ebin][cbin], hN1MCm22SUB2_narrow[ebin][cbin] );
            GraphToHist( gN1MCm18SUB2_narrow[ebin][cbin], hN1MCp18SUB2_narrow[ebin][cbin] );
            GraphToHist( gN1MCm14SUB2_narrow[ebin][cbin], hN1MCp14SUB2_narrow[ebin][cbin] );
            GraphToHist( gN1MCp22SUB2_narrow[ebin][cbin], hN1MCp22SUB2_narrow[ebin][cbin] );
            GraphToHist( gN1MCp18SUB2_narrow[ebin][cbin], hN1MCp18SUB2_narrow[ebin][cbin] );
            GraphToHist( gN1MCp14SUB2_narrow[ebin][cbin], hN1MCp14SUB2_narrow[ebin][cbin] );

            GraphToHist( gN1MCm22SUB3_narrow[ebin][cbin], hN1MCm22SUB3_narrow[ebin][cbin] );
            GraphToHist( gN1MCm18SUB3_narrow[ebin][cbin], hN1MCp18SUB3_narrow[ebin][cbin] );
            GraphToHist( gN1MCm14SUB3_narrow[ebin][cbin], hN1MCp14SUB3_narrow[ebin][cbin] );
            GraphToHist( gN1MCp22SUB3_narrow[ebin][cbin], hN1MCp22SUB3_narrow[ebin][cbin] );
            GraphToHist( gN1MCp18SUB3_narrow[ebin][cbin], hN1MCp18SUB3_narrow[ebin][cbin] );
            GraphToHist( gN1MCp14SUB3_narrow[ebin][cbin], hN1MCp14SUB3_narrow[ebin][cbin] );


            GraphToHist( gN1MCm22SUB2_wide[ebin][cbin], hN1MCm22SUB2_wide[ebin][cbin] );
            GraphToHist( gN1MCm18SUB2_wide[ebin][cbin], hN1MCp18SUB2_wide[ebin][cbin] );
            GraphToHist( gN1MCm14SUB2_wide[ebin][cbin], hN1MCp14SUB2_wide[ebin][cbin] );
            GraphToHist( gN1MCp22SUB2_wide[ebin][cbin], hN1MCp22SUB2_wide[ebin][cbin] );
            GraphToHist( gN1MCp18SUB2_wide[ebin][cbin], hN1MCp18SUB2_wide[ebin][cbin] );
            GraphToHist( gN1MCp14SUB2_wide[ebin][cbin], hN1MCp14SUB2_wide[ebin][cbin] );

            GraphToHist( gN1MCm22SUB3_wide[ebin][cbin], hN1MCm22SUB3_wide[ebin][cbin] );
            GraphToHist( gN1MCm18SUB3_wide[ebin][cbin], hN1MCp18SUB3_wide[ebin][cbin] );
            GraphToHist( gN1MCm14SUB3_wide[ebin][cbin], hN1MCp14SUB3_wide[ebin][cbin] );
            GraphToHist( gN1MCp22SUB3_wide[ebin][cbin], hN1MCp22SUB3_wide[ebin][cbin] );
            GraphToHist( gN1MCp18SUB3_wide[ebin][cbin], hN1MCp18SUB3_wide[ebin][cbin] );
            GraphToHist( gN1MCp14SUB3_wide[ebin][cbin], hN1MCp14SUB3_wide[ebin][cbin] );
        }
    }

}

void combineHists() {

    TH1::SetDefaultSumw2();

    tfin_nominal = new TFile("MH_hists.root","read");
    tfin_tight2 = new TFile("MH_tight2_hists.root","read");
    tfin_narrow = new TFile("MH_narrow_hists.root","read");
    tfin_wide = new TFile("MH_wide_hists.root","read");

    if (!fopen("results","r")) system("mkdir results");

    GetPtDists();

    GetEtaDists();

    tfin_nominal->Close();
    tfin_tight2->Close();
    tfin_narrow->Close();
    tfin_wide->Close();
}
