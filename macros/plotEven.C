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

TH1D * gN1MCm22SUB2_eta[ncentbins];
TH1D * gN1MCm18SUB2_eta[ncentbins];
TH1D * gN1MCm14SUB2_eta[ncentbins];
TH1D * gN1MCm10SUB2_eta[ncentbins];
TH1D * gN1MCm06SUB2_eta[ncentbins];
TH1D * gN1MCm02SUB2_eta[ncentbins];

TH1D * gN1MCp22SUB2_eta[ncentbins];
TH1D * gN1MCp18SUB2_eta[ncentbins];
TH1D * gN1MCp14SUB2_eta[ncentbins];
TH1D * gN1MCp10SUB2_eta[ncentbins];
TH1D * gN1MCp06SUB2_eta[ncentbins];
TH1D * gN1MCp02SUB2_eta[ncentbins];

TH1D * gN1MC22SUB2_eta[ncentbins];
TH1D * gN1MC18SUB2_eta[ncentbins];
TH1D * gN1MC14SUB2_eta[ncentbins];
TH1D * gN1MC10SUB2_eta[ncentbins];
TH1D * gN1MC06SUB2_eta[ncentbins];
TH1D * gN1MC02SUB2_eta[ncentbins];

TH1D * gN1MCm22SUB3_eta[ncentbins];
TH1D * gN1MCm18SUB3_eta[ncentbins];
TH1D * gN1MCm14SUB3_eta[ncentbins];
TH1D * gN1MCm10SUB3_eta[ncentbins];
TH1D * gN1MCm06SUB3_eta[ncentbins];
TH1D * gN1MCm02SUB3_eta[ncentbins];

TH1D * gN1MCp22SUB3_eta[ncentbins];
TH1D * gN1MCp18SUB3_eta[ncentbins];
TH1D * gN1MCp14SUB3_eta[ncentbins];
TH1D * gN1MCp10SUB3_eta[ncentbins];
TH1D * gN1MCp06SUB3_eta[ncentbins];
TH1D * gN1MCp02SUB3_eta[ncentbins];

TH1D * gN1MC22SUB3_eta[ncentbins];
TH1D * gN1MC18SUB3_eta[ncentbins];
TH1D * gN1MC14SUB3_eta[ncentbins];
TH1D * gN1MC10SUB3_eta[ncentbins];
TH1D * gN1MC06SUB3_eta[ncentbins];
TH1D * gN1MC02SUB3_eta[ncentbins];

void plotEven()
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

    for (int i = 0; i<NANALS; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString tag = Form("figures_MH/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
            // cout<<tag.Data()<<endl;
            ifstream fin(tag.Data());

            double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
            int neta = 0;
            while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
                // if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;

                if (i == 0) {
                    gN1MCm22SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm22SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 1) {
                    gN1MCm18SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm18SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 2) {
                    gN1MCm14SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm14SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 3) {
                    gN1MCm10SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm10SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 4) {
                    gN1MCm06SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm06SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 5) {
                    gN1MCm02SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm02SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 6) {
                    gN1MCp02SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp02SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 7) {
                    gN1MCp06SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp06SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 8) {
                    gN1MCp10SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp10SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 9) {
                    gN1MCp14SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp14SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 10) {
                    gN1MCp18SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp18SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 11) {
                    gN1MCp22SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp22SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);

                } else if (i == 12) {
                    gN1MCm22SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm22SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 13) {
                    gN1MCm18SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm18SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 14) {
                    gN1MCm14SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm14SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 15) {
                    gN1MCm10SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm10SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 16) {
                    gN1MCm06SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm06SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 17) {
                    gN1MCm02SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCm02SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 18) {
                    gN1MCp02SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp02SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 19) {
                    gN1MCp06SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp06SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 20) {
                    gN1MCp10SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp10SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 21) {
                    gN1MCp14SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp14SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 22) {
                    gN1MCp18SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp18SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 23) {
                    gN1MCp22SUB3_eta[cbin]->SetBinContent(neta+1, v1A);
                    gN1MCp22SUB3_eta[cbin]->SetBinError(neta+1, v1Ae);
                }  else {continue;}

                neta++;
            }

            gN1MC22SUB2_eta[cbin]->SetMarkerColor(kBlue);
            gN1MC18SUB2_eta[cbin]->SetMarkerColor(kRed);
            gN1MC14SUB2_eta[cbin]->SetMarkerColor(kGreen+2);
            gN1MC10SUB2_eta[cbin]->SetMarkerColor(kMagenta);
            gN1MC06SUB2_eta[cbin]->SetMarkerColor(kOrange+7);
            gN1MC02SUB2_eta[cbin]->SetMarkerColor(kBlack);

            gN1MC22SUB2_eta[cbin]->SetLineColor(kBlue);
            gN1MC18SUB2_eta[cbin]->SetLineColor(kRed);
            gN1MC14SUB2_eta[cbin]->SetLineColor(kGreen+2);
            gN1MC10SUB2_eta[cbin]->SetLineColor(kMagenta);
            gN1MC06SUB2_eta[cbin]->SetLineColor(kOrange+7);
            gN1MC02SUB2_eta[cbin]->SetLineColor(kBlack);

            gN1MC22SUB2_eta[cbin]->SetMarkerStyle(21);
            gN1MC18SUB2_eta[cbin]->SetMarkerStyle(20);
            gN1MC14SUB2_eta[cbin]->SetMarkerStyle(34);
            gN1MC10SUB2_eta[cbin]->SetMarkerStyle(33);
            gN1MC06SUB2_eta[cbin]->SetMarkerStyle(21);
            gN1MC02SUB2_eta[cbin]->SetMarkerStyle(34);

            gN1MC22SUB2_eta[cbin]->SetMarkerSize(1.2);
            gN1MC18SUB2_eta[cbin]->SetMarkerSize(1.3);
            gN1MC14SUB2_eta[cbin]->SetMarkerSize(1.7);
            gN1MC10SUB2_eta[cbin]->SetMarkerSize(1.8);
            gN1MC06SUB2_eta[cbin]->SetMarkerSize(1.2);
            gN1MC02SUB2_eta[cbin]->SetMarkerSize(1.7);
        }
    }

    double emin1 = -2.4;
    double emax1 = -0.0001;
    double emin2 = 0.001;
    double emax2 = 2.4;
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        for (int ebin = 1; ebin<=netabins; ebin++) {
            if (gN1MC22SUB2_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC22SUB2_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC22SUB2_eta[cbin]->SetBinContent(ebin, gN1MCp22SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC22SUB2_eta[cbin]->SetBinError(ebin, gN1MCp22SUB2_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC22SUB2_eta[cbin]->SetBinContent(ebin, gN1MCm22SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC22SUB2_eta[cbin]->SetBinError(ebin, gN1MCm22SUB2_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC18SUB2_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC18SUB2_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC18SUB2_eta[cbin]->SetBinContent(ebin, gN1MCp18SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC18SUB2_eta[cbin]->SetBinError(ebin, gN1MCp18SUB2_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC18SUB2_eta[cbin]->SetBinContent(ebin, gN1MCm18SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC18SUB2_eta[cbin]->SetBinError(ebin, gN1MCm18SUB2_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC14SUB2_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC14SUB2_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC14SUB2_eta[cbin]->SetBinContent(ebin, gN1MCp14SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC14SUB2_eta[cbin]->SetBinError(ebin, gN1MCp14SUB2_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC14SUB2_eta[cbin]->SetBinContent(ebin, gN1MCm14SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC14SUB2_eta[cbin]->SetBinError(ebin, gN1MCm14SUB2_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC10SUB2_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC10SUB2_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC10SUB2_eta[cbin]->SetBinContent(ebin, gN1MCp10SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC10SUB2_eta[cbin]->SetBinError(ebin, gN1MCp10SUB2_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC10SUB2_eta[cbin]->SetBinContent(ebin, gN1MCm10SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC10SUB2_eta[cbin]->SetBinError(ebin, gN1MCm10SUB2_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC06SUB2_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC06SUB2_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC06SUB2_eta[cbin]->SetBinContent(ebin, gN1MCp06SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC06SUB2_eta[cbin]->SetBinError(ebin, gN1MCp06SUB2_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC06SUB2_eta[cbin]->SetBinContent(ebin, gN1MCm06SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC06SUB2_eta[cbin]->SetBinError(ebin, gN1MCm06SUB2_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC02SUB2_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC02SUB2_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC02SUB2_eta[cbin]->SetBinContent(ebin, gN1MCp02SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC02SUB2_eta[cbin]->SetBinError(ebin, gN1MCp02SUB2_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC02SUB2_eta[cbin]->SetBinContent(ebin, gN1MCm02SUB2_eta[cbin]->GetBinContent(ebin));
                gN1MC02SUB2_eta[cbin]->SetBinError(ebin, gN1MCm02SUB2_eta[cbin]->GetBinError(ebin));
            }

            if (gN1MC22SUB3_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC22SUB3_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC22SUB3_eta[cbin]->SetBinContent(ebin, gN1MCp22SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC22SUB3_eta[cbin]->SetBinError(ebin, gN1MCp22SUB3_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC22SUB3_eta[cbin]->SetBinContent(ebin, gN1MCm22SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC22SUB3_eta[cbin]->SetBinError(ebin, gN1MCm22SUB3_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC18SUB3_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC18SUB3_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC18SUB3_eta[cbin]->SetBinContent(ebin, gN1MCp18SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC18SUB3_eta[cbin]->SetBinError(ebin, gN1MCp18SUB3_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC18SUB3_eta[cbin]->SetBinContent(ebin, gN1MCm18SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC18SUB3_eta[cbin]->SetBinError(ebin, gN1MCm18SUB3_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC14SUB3_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC14SUB3_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC14SUB3_eta[cbin]->SetBinContent(ebin, gN1MCp14SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC14SUB3_eta[cbin]->SetBinError(ebin, gN1MCp14SUB3_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC14SUB3_eta[cbin]->SetBinContent(ebin, gN1MCm14SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC14SUB3_eta[cbin]->SetBinError(ebin, gN1MCm14SUB3_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC10SUB3_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC10SUB3_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC10SUB3_eta[cbin]->SetBinContent(ebin, gN1MCp10SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC10SUB3_eta[cbin]->SetBinError(ebin, gN1MCp10SUB3_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC10SUB3_eta[cbin]->SetBinContent(ebin, gN1MCm10SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC10SUB3_eta[cbin]->SetBinError(ebin, gN1MCm10SUB3_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC06SUB3_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC06SUB3_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC06SUB3_eta[cbin]->SetBinContent(ebin, gN1MCp06SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC06SUB3_eta[cbin]->SetBinError(ebin, gN1MCp06SUB3_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC06SUB3_eta[cbin]->SetBinContent(ebin, gN1MCm06SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC06SUB3_eta[cbin]->SetBinError(ebin, gN1MCm06SUB3_eta[cbin]->GetBinError(ebin));
            }
            if (gN1MC02SUB3_eta[cbin]->GetBinCenter(ebin)>emin1 && gN1MC02SUB3_eta[cbin]->GetBinCenter(ebin)<emax1) {
                gN1MC02SUB3_eta[cbin]->SetBinContent(ebin, gN1MCp02SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC02SUB3_eta[cbin]->SetBinError(ebin, gN1MCp02SUB3_eta[cbin]->GetBinError(ebin));
            } else {
                gN1MC02SUB3_eta[cbin]->SetBinContent(ebin, gN1MCm02SUB3_eta[cbin]->GetBinContent(ebin));
                gN1MC02SUB3_eta[cbin]->SetBinError(ebin, gN1MCm02SUB3_eta[cbin]->GetBinError(ebin));
            }
        }
    }

    if (!fopen("figures_MH/v1even","r")) system("mkdir figures_MH/v1even");

    int setcent = 5;

    TCanvas * cEtaCompareEP = new TCanvas("cEtaCompareEP","cEtaCompareEP",650,600);
    TPad * padEtaCompareEP = (TPad *) cEtaCompareEP->cd();
    padEtaCompareEP->SetGrid();
    TH1D * hEtaCompareEP = new TH1D("hEtaCompareEP", "", 100, -2.5, 2.5);
    hEtaCompareEP->SetTitle("");
    hEtaCompareEP->SetStats(0);
    hEtaCompareEP->GetYaxis()->SetRangeUser(-0.06, 0.02);
    hEtaCompareEP->SetXTitle("#eta");
    hEtaCompareEP->SetYTitle("v_{1}^{even}");
    hEtaCompareEP->Draw();
    gN1MC22SUB2_eta[setcent]->Draw("same");
    gN1MC18SUB2_eta[setcent]->Draw("same");
    gN1MC14SUB2_eta[setcent]->Draw("same");
    gN1MC10SUB2_eta[setcent]->Draw("same");
    gN1MC06SUB2_eta[setcent]->Draw("same");
    gN1MC02SUB2_eta[setcent]->Draw("same");
    TLegend * legEtaCompareEP = new TLegend(0.19, 0.18, 0.46, 0.44);
    SetLegend(legEtaCompareEP,18);
    legEtaCompareEP->AddEntry(gN1MC22SUB2_eta[setcent],"N1MC22SUB2","p");
    legEtaCompareEP->AddEntry(gN1MC18SUB2_eta[setcent],"N1MC18SUB2","p");
    legEtaCompareEP->AddEntry(gN1MC14SUB2_eta[setcent],"N1MC14SUB2","p");
    legEtaCompareEP->AddEntry(gN1MC10SUB2_eta[setcent],"N1MC10SUB2","p");
    legEtaCompareEP->AddEntry(gN1MC06SUB2_eta[setcent],"N1MC06SUB2","p");
    legEtaCompareEP->AddEntry(gN1MC02SUB2_eta[setcent],"N1MC02SUB2","p");
    legEtaCompareEP->Draw();
    TPaveText * txEtaCompareEP = new TPaveText(0.20, 0.81, 0.47, 0.93, "NDC");
    SetTPaveTxt(txEtaCompareEP, 20);
    txEtaCompareEP->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txEtaCompareEP->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txEtaCompareEP->Draw();
    cEtaCompareEP->Print(Form("figures_MH/v1even/CompareEPrange_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");

}
