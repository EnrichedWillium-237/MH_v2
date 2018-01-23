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
static const int netabins = 14;
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0, -2.4,  0.0};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  0.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24", "-24_0", "0_24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCp22SUB2", "N1MCp18SUB2", "N1MCp14SUB2",
                  "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCp22SUB3", "N1MCp18SUB3", "N1MCp14SUB3",
                  "N1SUB2",      "N1ASUB2",     "N1BSUB2",     "N112SUB2",    "N112ASUB2",   "N112BSUB2",
                  "N1SUB3",      "N1ASUB3",     "N1BSUB3",     "N112SUB3",    "N112ASUB3",   "N112BSUB3"
};

TGraphErrors * gN1MCm22SUB2[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB2[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB2[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB2[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB2[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB2[netabins][ncentbins];
TGraphErrors * gN1SUB2[netabins][ncentbins];
TGraphErrors * gN1ASUB2[netabins][ncentbins];
TGraphErrors * gN1BSUB2[netabins][ncentbins];
TGraphErrors * gN112SUB2[netabins][ncentbins];
TGraphErrors * gN112ASUB2[netabins][ncentbins];
TGraphErrors * gN112BSUB2[netabins][ncentbins];

TGraphErrors * gN1MCm22SUB3[netabins][ncentbins];
TGraphErrors * gN1MCm18SUB3[netabins][ncentbins];
TGraphErrors * gN1MCm14SUB3[netabins][ncentbins];
TGraphErrors * gN1MCp22SUB3[netabins][ncentbins];
TGraphErrors * gN1MCp18SUB3[netabins][ncentbins];
TGraphErrors * gN1MCp14SUB3[netabins][ncentbins];
TGraphErrors * gN1SUB3[netabins][ncentbins];
TGraphErrors * gN1ASUB3[netabins][ncentbins];
TGraphErrors * gN1BSUB3[netabins][ncentbins];
TGraphErrors * gN112SUB3[netabins][ncentbins];
TGraphErrors * gN112ASUB3[netabins][ncentbins];
TGraphErrors * gN112BSUB3[netabins][ncentbins];

TH1D * hN1MCm22SUB2[netabins][ncentbins];
TH1D * hN1MCm18SUB2[netabins][ncentbins];
TH1D * hN1MCm14SUB2[netabins][ncentbins];
TH1D * hN1MCp22SUB2[netabins][ncentbins];
TH1D * hN1MCp18SUB2[netabins][ncentbins];
TH1D * hN1MCp14SUB2[netabins][ncentbins];
TH1D * hN1SUB2[netabins][ncentbins];
TH1D * hN1ASUB2[netabins][ncentbins];
TH1D * hN1BSUB2[netabins][ncentbins];
TH1D * hN112SUB2[netabins][ncentbins];
TH1D * hN112ASUB2[netabins][ncentbins];
TH1D * hN112BSUB2[netabins][ncentbins];

TH1D * hN1MCm22SUB3[netabins][ncentbins];
TH1D * hN1MCm18SUB3[netabins][ncentbins];
TH1D * hN1MCm14SUB3[netabins][ncentbins];
TH1D * hN1MCp22SUB3[netabins][ncentbins];
TH1D * hN1MCp18SUB3[netabins][ncentbins];
TH1D * hN1MCp14SUB3[netabins][ncentbins];
TH1D * hN1SUB3[netabins][ncentbins];
TH1D * hN1ASUB3[netabins][ncentbins];
TH1D * hN1BSUB3[netabins][ncentbins];
TH1D * hN112SUB3[netabins][ncentbins];
TH1D * hN112ASUB3[netabins][ncentbins];
TH1D * hN112BSUB3[netabins][ncentbins];

TFile * tfin;
TFile * tfout;

TDirectory * tdir;

bool isNominal = true;
bool isTight2 = false;
bool isWide = false;
bool isNarrow = false;

void GetPtDists( TString input = "MH" )
{
    isNominal = true;
    isTight2 = false;
    isWide = false;
    isNarrow = false;
    if (input.Contains("tight2")) isTight2 = true;
    if (input.Contains("wide")) isWide = true;
    if (input.Contains("narrow")) isNarrow = true;
    if (isTight2) {
        tfin = new TFile("MH_tight2_hists.root","read");
        tdir = (TDirectory *) tfout->mkdir("MH_tight2");
    } else if (isWide) {
        tfin = new TFile("MH_wide_hists.root","read");
        tdir = (TDirectory *) tfout->mkdir("MH_wide");
    } else if (isNarrow) {
        tfin = new TFile("MH_narrow_hists.root","read");
        tdir = (TDirectory *) tfout->mkdir("MH_narrow");
    } else {
        tfin = new TFile("MH_hists.root","read");
        tdir = (TDirectory *) tfout->mkdir("MH_nominal");
    }

    for (int i = 0; i<NANALS; i++) {
        for (int ebin = 0; ebin<netabins; ebin++) {
            for (int cbin = 0; cbin<ncentbins; cbin++) {
                TString mtag = Form("%s/%s/%d_%d",ANAL[i].Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);
                //cout<<mtag.Data()<<endl;
                if (i == 0) gN1MCm22SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 1) gN1MCm18SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 2) gN1MCm14SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 3) gN1MCp22SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 4) gN1MCp18SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 5) gN1MCp14SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));

                if (i == 6) gN1MCm22SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 7) gN1MCm18SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 8) gN1MCm14SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 9) gN1MCp22SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 10) gN1MCp18SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 11) gN1MCp14SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));

                if (i == 12) gN1SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 13) gN1ASUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 14) gN1BSUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 15) gN112SUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 16) gN112ASUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 17) gN112BSUB2[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));

                if (i == 18) gN1SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 19) gN1ASUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 20) gN1BSUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 21) gN112SUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 22) gN112ASUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
                if (i == 23) gN112BSUB3[ebin][cbin] = (TGraphErrors *) tfin->Get(Form("%s/gA",mtag.Data()));
            }
        }
    }

    if (isTight2) cout << "Accessing MH_tight2_hists.root... " << endl;
    else if (isWide) cout << "Accessing MH_wide_hists.root... " << endl;
    else if (isNarrow) cout << "Accessing MH_narrow_hists.root... " << endl;
    else cout << "Accessing MH_hists.root... " << endl;

    for (int ebin = 0; ebin<netabins; ebin++) {
        TDirectory * tdebin = (TDirectory *) tdir->mkdir(Form("eta_%s",etags[ebin].Data()));
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * tdcbin = (TDirectory *) tdebin->mkdir(Form("%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
            tdcbin->cd();

            hN1MCm22SUB2[ebin][cbin] = new TH1D(Form("N1MCm22SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            hN1MCm22SUB2[ebin][cbin]->SetStats(0);
            hN1MCm22SUB2[ebin][cbin]->SetXTitle("p_{T} (GeV/c)");
            hN1MCm22SUB2[ebin][cbin]->SetYTitle("v_{1}");

            hN1MCm18SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm18SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm14SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp22SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp18SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp14SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            hN1MCm22SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm22SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm18SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm18SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCm14SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCm14SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp22SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp22SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp18SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp18SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1MCp14SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1MCp14SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            hN1SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1ASUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1ASUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1BSUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1BSUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN112SUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N112SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN112ASUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N112ASUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN112BSUB2[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N112BSUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            hN1SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1ASUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1ASUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN1BSUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N1BSUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN112SUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N112SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN112ASUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N112ASUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            hN112BSUB3[ebin][cbin] = (TH1D *) hN1MCm22SUB2[ebin][cbin]->Clone(Form("N112BSUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            // convert TGraphErrors to TH1Ds
            GraphToHist( gN1MCm22SUB2[ebin][cbin], hN1MCm22SUB2[ebin][cbin] );
            GraphToHist( gN1MCm18SUB2[ebin][cbin], hN1MCm18SUB2[ebin][cbin] );
            GraphToHist( gN1MCm14SUB2[ebin][cbin], hN1MCm14SUB2[ebin][cbin] );
            GraphToHist( gN1MCp22SUB2[ebin][cbin], hN1MCp22SUB2[ebin][cbin] );
            GraphToHist( gN1MCp18SUB2[ebin][cbin], hN1MCp18SUB2[ebin][cbin] );
            GraphToHist( gN1MCp14SUB2[ebin][cbin], hN1MCp14SUB2[ebin][cbin] );

            GraphToHist( gN1MCm22SUB3[ebin][cbin], hN1MCm22SUB3[ebin][cbin] );
            GraphToHist( gN1MCm18SUB3[ebin][cbin], hN1MCm18SUB3[ebin][cbin] );
            GraphToHist( gN1MCm14SUB3[ebin][cbin], hN1MCm14SUB3[ebin][cbin] );
            GraphToHist( gN1MCp22SUB3[ebin][cbin], hN1MCp22SUB3[ebin][cbin] );
            GraphToHist( gN1MCp18SUB3[ebin][cbin], hN1MCp18SUB3[ebin][cbin] );
            GraphToHist( gN1MCp14SUB3[ebin][cbin], hN1MCp14SUB3[ebin][cbin] );

            GraphToHist( gN1SUB2[ebin][cbin], hN1SUB2[ebin][cbin] );
            GraphToHist( gN1ASUB2[ebin][cbin], hN1ASUB2[ebin][cbin] );
            GraphToHist( gN1BSUB2[ebin][cbin], hN1BSUB2[ebin][cbin] );
            GraphToHist( gN112SUB2[ebin][cbin], hN112SUB2[ebin][cbin] );
            GraphToHist( gN112ASUB2[ebin][cbin], hN112ASUB2[ebin][cbin] );
            GraphToHist( gN112BSUB2[ebin][cbin], hN112BSUB2[ebin][cbin] );

            GraphToHist( gN1SUB3[ebin][cbin], hN1SUB3[ebin][cbin] );
            GraphToHist( gN1ASUB3[ebin][cbin], hN1ASUB3[ebin][cbin] );
            GraphToHist( gN1BSUB3[ebin][cbin], hN1BSUB3[ebin][cbin] );
            GraphToHist( gN112SUB3[ebin][cbin], hN112SUB3[ebin][cbin] );
            GraphToHist( gN112ASUB3[ebin][cbin], hN112ASUB3[ebin][cbin] );
            GraphToHist( gN112BSUB3[ebin][cbin], hN112BSUB3[ebin][cbin] );

            // write histograms to output file
            hN1MCm22SUB2[ebin][cbin]->Write();
            hN1MCm18SUB2[ebin][cbin]->Write();
            hN1MCm14SUB2[ebin][cbin]->Write();
            hN1MCp22SUB2[ebin][cbin]->Write();
            hN1MCp18SUB2[ebin][cbin]->Write();
            hN1MCp14SUB2[ebin][cbin]->Write();

            hN1MCm22SUB3[ebin][cbin]->Write();
            hN1MCm18SUB3[ebin][cbin]->Write();
            hN1MCm14SUB3[ebin][cbin]->Write();
            hN1MCp22SUB3[ebin][cbin]->Write();
            hN1MCp18SUB3[ebin][cbin]->Write();
            hN1MCp14SUB3[ebin][cbin]->Write();

            hN1SUB2[ebin][cbin]->Write();
            hN1ASUB2[ebin][cbin]->Write();
            hN1BSUB2[ebin][cbin]->Write();
            hN112SUB2[ebin][cbin]->Write();
            hN112ASUB2[ebin][cbin]->Write();
            hN112BSUB2[ebin][cbin]->Write();

            hN1SUB3[ebin][cbin]->Write();
            hN1ASUB3[ebin][cbin]->Write();
            hN1BSUB3[ebin][cbin]->Write();
            hN112SUB3[ebin][cbin]->Write();
            hN112ASUB3[ebin][cbin]->Write();
            hN112BSUB3[ebin][cbin]->Write();

            // clean up histograms
            hN1MCm22SUB2[ebin][cbin]->Delete();
            hN1MCm18SUB2[ebin][cbin]->Delete();
            hN1MCm14SUB2[ebin][cbin]->Delete();
            hN1MCp22SUB2[ebin][cbin]->Delete();
            hN1MCp18SUB2[ebin][cbin]->Delete();
            hN1MCp14SUB2[ebin][cbin]->Delete();

            hN1MCm22SUB3[ebin][cbin]->Delete();
            hN1MCm18SUB3[ebin][cbin]->Delete();
            hN1MCm14SUB3[ebin][cbin]->Delete();
            hN1MCp22SUB3[ebin][cbin]->Delete();
            hN1MCp18SUB3[ebin][cbin]->Delete();
            hN1MCp14SUB3[ebin][cbin]->Delete();

            hN1SUB2[ebin][cbin]->Delete();
            hN1ASUB2[ebin][cbin]->Delete();
            hN1BSUB2[ebin][cbin]->Delete();
            hN112SUB2[ebin][cbin]->Delete();
            hN112ASUB2[ebin][cbin]->Delete();
            hN112BSUB2[ebin][cbin]->Delete();

            hN1SUB3[ebin][cbin]->Delete();
            hN1ASUB3[ebin][cbin]->Delete();
            hN1BSUB3[ebin][cbin]->Delete();
            hN112SUB3[ebin][cbin]->Delete();
            hN112ASUB3[ebin][cbin]->Delete();
            hN112BSUB3[ebin][cbin]->Delete();
        }
    }

    tfin->Close();

}

void combinePt() {

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
