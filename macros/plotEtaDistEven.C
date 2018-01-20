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

void plotEtaDistEven()
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

    if (!fopen("figures_MH/EtaDists","r")) system("mkdir figures_MH/EtaDists");

    int setcent = 5;
/*
    TCanvas * cEtaCompareEP = new TCanvas("cEtaCompareEP","cEtaCompareEP",650,600);
    TPad * padEtaCompareEP = (TPad *) cEtaCompareEP->cd();
    padEtaCompareEP->SetGrid();
    TH1D * hEtaCompareEP = new TH1D("hEtaCompareEP", "", 100, -2.5, 2.5);
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
    cEtaCompareEP->Print(Form("figures_MH/EtaDists/CompareEPrange_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cEtaMC22_MC18 = new TCanvas("cEtaMC22_MC18","cEtaMC22_MC18",650,600);
    TPad * padEtaMC22_MC18 = (TPad *) cEtaMC22_MC18->cd();
    padEtaMC22_MC18->SetGrid();
    TH1D * hEtaMC22_MC18 = new TH1D("hEtaMC22_MC18", "hEtaMC22_MC18", 100, -2.5, 2.5);
    hEtaMC22_MC18->SetStats(0);
    hEtaMC22_MC18->GetYaxis()->SetRangeUser(-0.02, 0.01);
    hEtaMC22_MC18->SetXTitle("#eta");
    hEtaMC22_MC18->SetYTitle("v_{1}^{even}");
    hEtaMC22_MC18->Draw();
    gN1MC22SUB2_eta[setcent]->Draw("same");
    gN1MC18SUB2_eta[setcent]->SetMarkerStyle(25);
    gN1MC18SUB2_eta[setcent]->SetMarkerSize(1.2);
    gN1MC18SUB2_eta[setcent]->Draw("same");
    TLegend * legEtaMC22_MC18 = new TLegend(0.65, 0.81, 0.93, 0.93);
    SetLegend(legEtaMC22_MC18, 18);
    legEtaMC22_MC18->AddEntry(gN1MC22SUB2_eta[setcent],"N1MC22SUB2","p");
    legEtaMC22_MC18->AddEntry(gN1MC18SUB2_eta[setcent],"N1MC18SUB2","p");
    legEtaMC22_MC18->Draw();
    TPaveText * txEtaMC22_MC18 = new TPaveText(0.20, 0.81, 0.47, 0.93, "NDC");
    SetTPaveTxt(txEtaMC22_MC18, 20);
    txEtaMC22_MC18->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txEtaMC22_MC18->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txEtaMC22_MC18->Draw();
    cEtaMC22_MC18->Print(Form("figures_MH/EtaDists/CompareMC22_MC18_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cEtaMC18_MC14 = new TCanvas("cEtaMC18_MC14","cEtaMC18_MC14",650,600);
    TPad * padEtaMC18_MC14 = (TPad *) cEtaMC18_MC14->cd();
    padEtaMC18_MC14->SetGrid();
    TH1D * hEtaMC18_MC14 = new TH1D("hEtaMC18_MC14", "hEtaMC18_MC14", 100, -2.5, 2.5);
    hEtaMC18_MC14->SetStats(0);
    hEtaMC18_MC14->GetYaxis()->SetRangeUser(-0.02, 0.01);
    hEtaMC18_MC14->SetXTitle("#eta");
    hEtaMC18_MC14->SetYTitle("v_{1}^{even}");
    hEtaMC18_MC14->Draw();
    gN1MC18SUB2_eta[setcent]->Draw("same");
    gN1MC14SUB2_eta[setcent]->SetMarkerStyle(28);
    gN1MC14SUB2_eta[setcent]->SetMarkerSize(1.7);
    gN1MC14SUB2_eta[setcent]->Draw("same");
    TLegend * legEtaMC18_MC14 = new TLegend(0.65, 0.81, 0.93, 0.93);
    SetLegend(legEtaMC18_MC14, 18);
    legEtaMC18_MC14->AddEntry(gN1MC18SUB2_eta[setcent],"N1MC18SUB2","p");
    legEtaMC18_MC14->AddEntry(gN1MC14SUB2_eta[setcent],"N1MC14SUB2","p");
    legEtaMC18_MC14->Draw();
    TPaveText * txEtaMC18_MC14 = new TPaveText(0.20, 0.81, 0.47, 0.93, "NDC");
    SetTPaveTxt(txEtaMC18_MC14, 20);
    txEtaMC18_MC14->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txEtaMC18_MC14->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txEtaMC18_MC14->Draw();
    cEtaMC18_MC14->Print(Form("figures_MH/EtaDists/CompareMC18_MC14_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cMCCent22eta = new TCanvas("cMCCent22eta","cMCCent22eta",650,600);
    TPad * padMCCent22eta = (TPad *) cMCCent22eta->cd();
    padMCCent22eta->SetGrid();
    TH1D * hMCCent22eta = new TH1D("hMCCent22eta", "", 100, -2.5, 2.5);
    hMCCent22eta->SetTitle("");
    hMCCent22eta->SetStats(0);
    hMCCent22eta->GetYaxis()->SetRangeUser(-0.015, 0.015);
    hMCCent22eta->SetXTitle("#eta");
    hMCCent22eta->SetYTitle("v_{1}^{even}");
    hMCCent22eta->Draw();
    TH1D * gN1MC22SUB2_tmp_0 = (TH1D *) gN1MC22SUB2_eta[0]->Clone("gN1MC22SUB2_eta_tmp_0");
    TH1D * gN1MC22SUB2_tmp_7 = (TH1D *) gN1MC22SUB2_eta[7]->Clone("gN1MC22SUB2_eta_tmp_7");
    TH1D * gN1MC22SUB2_tmp_8 = (TH1D *) gN1MC22SUB2_eta[8]->Clone("gN1MC22SUB2_eta_tmp_8");
    TH1D * gN1MC22SUB2_tmp_9 = (TH1D *) gN1MC22SUB2_eta[9]->Clone("gN1MC22SUB2_eta_tmp_9");
    TH1D * gN1MC22SUB2_tmp_10 = (TH1D *) gN1MC22SUB2_eta[10]->Clone("gN1MC22SUB2_eta_tmp_10");
    gN1MC22SUB2_tmp_0->SetMarkerColor(kBlue);
    gN1MC22SUB2_tmp_7->SetMarkerColor(kRed);
    gN1MC22SUB2_tmp_8->SetMarkerColor(kGreen+2);
    gN1MC22SUB2_tmp_9->SetMarkerColor(kMagenta);
    gN1MC22SUB2_tmp_10->SetMarkerColor(kOrange+7);
    gN1MC22SUB2_tmp_0->SetLineColor(kBlue);
    gN1MC22SUB2_tmp_7->SetLineColor(kRed);
    gN1MC22SUB2_tmp_8->SetLineColor(kGreen+2);
    gN1MC22SUB2_tmp_9->SetLineColor(kMagenta);
    gN1MC22SUB2_tmp_10->SetLineColor(kOrange+7);
    gN1MC22SUB2_tmp_0->SetMarkerStyle(21);
    gN1MC22SUB2_tmp_7->SetMarkerStyle(27);
    gN1MC22SUB2_tmp_8->SetMarkerStyle(34);
    gN1MC22SUB2_tmp_9->SetMarkerStyle(24);
    gN1MC22SUB2_tmp_10->SetMarkerStyle(33);
    gN1MC22SUB2_tmp_0->SetMarkerSize(1.2);
    gN1MC22SUB2_tmp_7->SetMarkerSize(1.8);
    gN1MC22SUB2_tmp_8->SetMarkerSize(1.7);
    gN1MC22SUB2_tmp_9->SetMarkerSize(1.3);
    gN1MC22SUB2_tmp_10->SetMarkerSize(1.8);
    gN1MC22SUB2_tmp_0->Draw("same");
    gN1MC22SUB2_tmp_7->Draw("same");
    gN1MC22SUB2_tmp_8->Draw("same");
    gN1MC22SUB2_tmp_9->Draw("same");
    //gN1MC22SUB2_tmp_10->Draw("same");
    TLegend * legMCCent22eta = new TLegend(0.74, 0.74, 0.93, 0.93);
    SetLegend(legMCCent22eta, 18);
    legMCCent22eta->AddEntry(gN1MC22SUB2_tmp_0,Form(" %d-%d%%",cminCENT[0],cmaxCENT[0]),"p");
    legMCCent22eta->AddEntry(gN1MC22SUB2_tmp_7,Form(" %d-%d%%",cminCENT[7],cmaxCENT[7]),"p");
    legMCCent22eta->AddEntry(gN1MC22SUB2_tmp_8,Form(" %d-%d%%",cminCENT[8],cmaxCENT[8]),"p");
    legMCCent22eta->AddEntry(gN1MC22SUB2_tmp_9,Form(" %d-%d%%",cminCENT[9],cmaxCENT[9]),"p");
    //legMCCent22eta->AddEntry(gN1MC22SUB2_tmp_10,Form(" %d-%d%%",cminCENT[10],cmaxCENT[10]),"p");
    legMCCent22eta->Draw();
    TPaveText * txMCCent22eta = new TPaveText(0.2, 0.81, 0.65, 0.92, "NDC");
    SetTPaveTxt(txMCCent22eta, 20);
    txMCCent22eta->AddText("#pm Event plane 2.0 < |#eta| < 2.4");
    txMCCent22eta->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCCent22eta->Draw();
    cMCCent22eta->Print(Form("figures_MH/EtaDists/CentScanMC22_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cMCCent18eta = new TCanvas("cMCCent18eta","cMCCent18eta",650,600);
    TPad * padMCCent18eta = (TPad *) cMCCent18eta->cd();
    padMCCent18eta->SetGrid();
    TH1D * hMCCent18eta = new TH1D("hMCCent18eta", "", 100, -2.5, 2.5);
    hMCCent18eta->SetTitle("");
    hMCCent18eta->SetStats(0);
    hMCCent18eta->GetYaxis()->SetRangeUser(-0.015, 0.015);
    hMCCent18eta->SetXTitle("#eta");
    hMCCent18eta->SetYTitle("v_{1}^{even}");
    hMCCent18eta->Draw();
    TH1D * gN1MC18SUB2_tmp_0 = (TH1D *) gN1MC18SUB2_eta[0]->Clone("gN1MC18SUB2_eta_tmp_0");
    TH1D * gN1MC18SUB2_tmp_7 = (TH1D *) gN1MC18SUB2_eta[7]->Clone("gN1MC18SUB2_eta_tmp_7");
    TH1D * gN1MC18SUB2_tmp_8 = (TH1D *) gN1MC18SUB2_eta[8]->Clone("gN1MC18SUB2_eta_tmp_8");
    TH1D * gN1MC18SUB2_tmp_9 = (TH1D *) gN1MC18SUB2_eta[9]->Clone("gN1MC18SUB2_eta_tmp_9");
    TH1D * gN1MC18SUB2_tmp_10 = (TH1D *) gN1MC18SUB2_eta[10]->Clone("gN1MC18SUB2_eta_tmp_10");
    gN1MC18SUB2_tmp_0->SetMarkerColor(kBlue);
    gN1MC18SUB2_tmp_7->SetMarkerColor(kRed);
    gN1MC18SUB2_tmp_8->SetMarkerColor(kGreen+2);
    gN1MC18SUB2_tmp_9->SetMarkerColor(kMagenta);
    gN1MC18SUB2_tmp_10->SetMarkerColor(kOrange+7);
    gN1MC18SUB2_tmp_0->SetLineColor(kBlue);
    gN1MC18SUB2_tmp_7->SetLineColor(kRed);
    gN1MC18SUB2_tmp_8->SetLineColor(kGreen+2);
    gN1MC18SUB2_tmp_9->SetLineColor(kMagenta);
    gN1MC18SUB2_tmp_10->SetLineColor(kOrange+7);
    gN1MC18SUB2_tmp_0->SetMarkerStyle(21);
    gN1MC18SUB2_tmp_7->SetMarkerStyle(27);
    gN1MC18SUB2_tmp_8->SetMarkerStyle(34);
    gN1MC18SUB2_tmp_9->SetMarkerStyle(24);
    gN1MC18SUB2_tmp_10->SetMarkerStyle(33);
    gN1MC18SUB2_tmp_0->SetMarkerSize(1.2);
    gN1MC18SUB2_tmp_7->SetMarkerSize(1.8);
    gN1MC18SUB2_tmp_8->SetMarkerSize(1.7);
    gN1MC18SUB2_tmp_9->SetMarkerSize(1.3);
    gN1MC18SUB2_tmp_10->SetMarkerSize(1.8);
    gN1MC18SUB2_tmp_0->Draw("same");
    gN1MC18SUB2_tmp_7->Draw("same");
    gN1MC18SUB2_tmp_8->Draw("same");
    gN1MC18SUB2_tmp_9->Draw("same");
    //gN1MC18SUB2_tmp_10->Draw("same");
    TLegend * legMCCent18eta = new TLegend(0.74, 0.74, 0.93, 0.93);
    SetLegend(legMCCent18eta, 18);
    legMCCent18eta->AddEntry(gN1MC18SUB2_tmp_0,Form(" %d-%d%%",cminCENT[0],cmaxCENT[0]),"p");
    legMCCent18eta->AddEntry(gN1MC18SUB2_tmp_7,Form(" %d-%d%%",cminCENT[7],cmaxCENT[7]),"p");
    legMCCent18eta->AddEntry(gN1MC18SUB2_tmp_8,Form(" %d-%d%%",cminCENT[8],cmaxCENT[8]),"p");
    legMCCent18eta->AddEntry(gN1MC18SUB2_tmp_9,Form(" %d-%d%%",cminCENT[9],cmaxCENT[9]),"p");
    //legMCCent18eta->AddEntry(gN1MC18SUB2_tmp_10,Form(" %d-%d%%",cminCENT[10],cmaxCENT[10]),"p");
    legMCCent18eta->Draw();
    TPaveText * txMCCent18eta = new TPaveText(0.2, 0.81, 0.65, 0.92, "NDC");
    SetTPaveTxt(txMCCent18eta, 20);
    txMCCent18eta->AddText("#pm Event plane 1.6 < |#eta| < 2.0");
    txMCCent18eta->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCCent18eta->Draw();
    cMCCent18eta->Print(Form("figures_MH/EtaDists/CentScanMC18_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cEtaCompareSUB_MC22 = new TCanvas("cEtaCompareSUB_MC22","cEtaCompareSUB_MC22",650,600);
    TPad * padEtaCompareSUB_MC22 = (TPad *) cEtaCompareSUB_MC22->cd();
    padEtaCompareSUB_MC22->SetGrid();
    TH1D * hEtaCompareSUB_MC22 = new TH1D("hEtaCompareSUB_MC22", "hEtaCompareSUB_MC22", 100, -2.5, 2.5);
    hEtaCompareSUB_MC22->SetStats(0);
    hEtaCompareSUB_MC22->GetYaxis()->SetRangeUser(-0.02, 0.01);
    hEtaCompareSUB_MC22->SetXTitle("#eta");
    hEtaCompareSUB_MC22->SetYTitle("v_{1}^{even}");
    hEtaCompareSUB_MC22->Draw();
    gN1MC22SUB2_eta[setcent]->Draw("same");
    gN1MC22SUB3_eta[setcent]->SetMarkerStyle(25);
    gN1MC22SUB3_eta[setcent]->SetMarkerSize(1.2);
    gN1MC22SUB3_eta[setcent]->Draw("same");
    TLegend * legEtaCompareSUB_MC22 = new TLegend(0.65, 0.81, 0.93, 0.93);
    SetLegend(legEtaCompareSUB_MC22, 18);
    legEtaCompareSUB_MC22->AddEntry(gN1MC22SUB2_eta[setcent],"N1MC22SUB2","p");
    legEtaCompareSUB_MC22->AddEntry(gN1MC22SUB3_eta[setcent],"N1MC22SUB3","p");
    legEtaCompareSUB_MC22->Draw();
    TPaveText * txEtaCompareSUB_MC22 = new TPaveText(0.20, 0.81, 0.47, 0.93, "NDC");
    SetTPaveTxt(txEtaCompareSUB_MC22, 20);
    txEtaCompareSUB_MC22->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txEtaCompareSUB_MC22->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txEtaCompareSUB_MC22->Draw();
    cEtaCompareSUB_MC22->Print(Form("figures_MH/EtaDists/CompareSUB_MC22_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");
*/

    TCanvas * cEtaCompareSUB_MC18 = new TCanvas("cEtaCompareSUB_MC18","cEtaCompareSUB_MC18",650,600);
    TPad * padEtaCompareSUB_MC18 = (TPad *) cEtaCompareSUB_MC18->cd();
    padEtaCompareSUB_MC18->SetGrid();
    TH1D * hEtaCompareSUB_MC18 = new TH1D("hEtaCompareSUB_MC18", "", 100, -2.5, 2.5);
    hEtaCompareSUB_MC18->SetStats(0);
    hEtaCompareSUB_MC18->GetYaxis()->SetRangeUser(-0.02, 0.01);
    hEtaCompareSUB_MC18->SetXTitle("#eta");
    hEtaCompareSUB_MC18->SetYTitle("v_{1}^{even}");
    hEtaCompareSUB_MC18->Draw();
    gN1MC18SUB2_eta[setcent]->Draw("same");
    gN1MC18SUB3_eta[setcent]->SetMarkerStyle(25);
    gN1MC18SUB3_eta[setcent]->SetMarkerSize(1.2);
    gN1MC18SUB3_eta[setcent]->Draw("same");
    TLegend * legEtaCompareSUB_MC18 = new TLegend(0.65, 0.81, 0.93, 0.93);
    SetLegend(legEtaCompareSUB_MC18, 18);
    legEtaCompareSUB_MC18->AddEntry(gN1MC18SUB2_eta[setcent],"N1MC18SUB2","p");
    legEtaCompareSUB_MC18->AddEntry(gN1MC18SUB3_eta[setcent],"N1MC18SUB3","p");
    legEtaCompareSUB_MC18->Draw();
    TPaveText * txEtaCompareSUB_MC18 = new TPaveText(0.20, 0.81, 0.47, 0.93, "NDC");
    SetTPaveTxt(txEtaCompareSUB_MC18, 20);
    txEtaCompareSUB_MC18->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txEtaCompareSUB_MC18->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txEtaCompareSUB_MC18->Draw();
    cEtaCompareSUB_MC18->Print(Form("figures_MH/EtaDists/CompareSUB_MC18_eta_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cEtapmMC22 = new TCanvas("cEtapmMC22","cEtapmMC22",800,600);
    TPad * padEtapmMC22 = (TPad *) cEtapmMC22->cd();
    padEtapmMC22->SetGrid();
    TH1D * hEtapmMC22 = new TH1D("hEtapmMC22", "", 100, -2.5, 2.5);
    hEtapmMC22->SetStats(0);
    hEtapmMC22->GetYaxis()->SetRangeUser(-0.1, 0.2);
    hEtapmMC22->SetXTitle("#eta");
    hEtapmMC22->SetYTitle("v_{1}^{even}");
    hEtapmMC22->Draw("same");
    gN1MCm22SUB2_eta[setcent]->SetMarkerColor(kRed);
    gN1MCm22SUB2_eta[setcent]->SetLineColor(kRed);
    gN1MCm22SUB2_eta[setcent]->SetMarkerStyle(24);
    gN1MCm22SUB2_eta[setcent]->SetMarkerSize(1.2);
    gN1MCm22SUB2_eta[setcent]->Draw("same");
    gN1MCp22SUB2_eta[setcent]->SetMarkerColor(kBlue);
    gN1MCp22SUB2_eta[setcent]->SetLineColor(kBlue);
    gN1MCp22SUB2_eta[setcent]->SetMarkerStyle(24);
    gN1MCp22SUB2_eta[setcent]->SetMarkerSize(1.2);
    gN1MCp22SUB2_eta[setcent]->Draw("same");
    TLegend * legEtapmMC22 = new TLegend(0.2, 0.2, 0.4, 0.4);
    SetLegend(legEtapmMC22, 18);
    legEtapmMC22->AddEntry(gN1MCm22SUB2_eta[setcent],"#Psi_{1}^{trk} (-2.4 < #eta < -2.0)","p");
    legEtapmMC22->AddEntry(gN1MCp22SUB2_eta[setcent],"#Psi_{1}^{trk} (2.0 < #eta < 2.4)","p");
    legEtapmMC22->Draw();
    TPaveText * txEtapmMC22 = new TPaveText(0.4, 0.4, 0.6, 0.6, "NDC");
    SetTPaveTxt(xEtapmMC22, 20);
    txEtapmMC22->AddEntry("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txEtapmMC22->AddEntry(Form("(%d-%d%%)",cminCENT[setcent],cmaxCENT[setcent]));
    txEtapmMC22->Draw();
    cEtapmMC22->Print(Form("EtapmMC22_%d_%d",cminCENT[setcent],cmaxCENT[setcent]));


    TCanvas * cEtapmMC18 = new TCanvas("cEtapmMC18","cEtapmMC18",800,600);
    TPad * padEtapmMC18 = (TPad *) cEtapmMC18->cd();
    padEtapmMC18->SetGrid();
    TH1D * hEtapmMC18 = new TH1D("hEtapmMC18", "", 100, -2.5, 2.5);
    hEtapmMC18->SetStats(0);
    hEtapmMC18->GetYaxis()->SetRangeUser(-0.1, 0.2);
    hEtapmMC18->SetXTitle("#eta");
    hEtapmMC18->SetYTitle("v_{1}^{even}");
    hEtapmMC18->Draw("same");
    gN1MCm18SUB2_eta[setcent]->SetMarkerColor(kRed);
    gN1MCm18SUB2_eta[setcent]->SetLineColor(kRed);
    gN1MCm18SUB2_eta[setcent]->SetMarkerStyle(24);
    gN1MCm18SUB2_eta[setcent]->SetMarkerSize(1.2);
    gN1MCm18SUB2_eta[setcent]->Draw("same");
    gN1MCp18SUB2_eta[setcent]->SetMarkerColor(kBlue);
    gN1MCp18SUB2_eta[setcent]->SetLineColor(kBlue);
    gN1MCp18SUB2_eta[setcent]->SetMarkerStyle(24);
    gN1MCp18SUB2_eta[setcent]->SetMarkerSize(1.2);
    gN1MCp18SUB2_eta[setcent]->Draw("same");
    TLegend * legEtapmMC18 = new TLegend(0.2, 0.2, 0.4, 0.4);
    SetLegend(legEtapmMC18, 18);
    legEtapmMC18->AddEntry(gN1MCm18SUB2_eta[setcent],"#Psi_{1}^{trk} (-2.0 < #eta < -1.6)","p");
    legEtapmMC18->AddEntry(gN1MCp18SUB2_eta[setcent],"#Psi_{1}^{trk} (1.6 < #eta < 2.0)","p");
    legEtapmMC18->Draw();
    TPaveText * txEtapmMC18 = new TPaveText(0.4, 0.4, 0.6, 0.6, "NDC");
    SetTPaveTxt(xEtapmMC18, 20);
    txEtapmMC18->AddEntry("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txEtapmMC18->AddEntry(Form("(%d-%d%%)",cminCENT[setcent],cmaxCENT[setcent]));
    txEtapmMC18->Draw();
    cEtapmMC18->Print(Form("EtapmMC18_%d_%d",cminCENT[setcent],cmaxCENT[setcent]));


}
