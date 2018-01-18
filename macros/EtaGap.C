# include "TArrow.h"
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

static const int cbinsCENT = 13;
static const int cminCENT[] = {0,  5, 10, 15, 20, 25, 30, 35, 40, 50, 60,  0, 20,  60};
static const int cmaxCENT[] = {5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 20, 60, 100};
static const int ndeltaEta = 11;
static const double deltaEta[] = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
static int NANALS = 24;
string ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCm10SUB2", "N1MCm06SUB2", "N1MCm02SUB2",
                 "N1MCp02SUB2", "N1MCp06SUB2", "N1MCp10SUB2", "N1MCp14SUB2", "N1MCp18SUB2", "N1MCp22SUB2",
                 "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCm10SUB3", "N1MCm06SUB3", "N1MCm02SUB3",
                 "N1MCp02SUB3", "N1MCp06SUB3", "N1MCp10SUB3", "N1MCp14SUB3", "N1MCp18SUB3", "N1MCp22SUB3"
};

TGraphErrors * gN1MCm22SUB2[cbinsCENT];
TGraphErrors * gN1MCm18SUB2[cbinsCENT];
TGraphErrors * gN1MCm14SUB2[cbinsCENT];
TGraphErrors * gN1MCm10SUB2[cbinsCENT];
TGraphErrors * gN1MCm06SUB2[cbinsCENT];
TGraphErrors * gN1MCm02SUB2[cbinsCENT];
TGraphErrors * gN1MCp02SUB2[cbinsCENT];
TGraphErrors * gN1MCp06SUB2[cbinsCENT];
TGraphErrors * gN1MCp10SUB2[cbinsCENT];
TGraphErrors * gN1MCp14SUB2[cbinsCENT];
TGraphErrors * gN1MCp18SUB2[cbinsCENT];
TGraphErrors * gN1MCp22SUB2[cbinsCENT];

void EtaGap()
{

    TH1::SetDefaultSumw2();

    for (int i = 0; i<NANALS; i++) {

        for (int cbin = 0; cbin<cbinsCENT; cbin++) {
            TString tag = Form("figures_MH/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].data(),ANAL[i].data(),cminCENT[cbin],cmaxCENT[cbin]);
            // cout<<tag.Data()<<endl;
            ifstream fin(tag.Data());

            double etaval[20] = {0};
            double v1val[20] = {0};
            double v1vale[20] = {0};
            double v1Aval[20] = {0};
            double v1Avale[20] = {0};
            double v1Bval[20] = {0};
            double v1Bvale[20] = {0};
            double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
            int neta = 0;
            while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
                //if (i == 0 && cbin == 0) cout<<"eta = "<<eta<<"\tv1 = "<<v1<<"\tv1e = "<<v1e<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<"\tv1B = "<<v1B<<"\tv1Be = "<<v1Be<<endl;
                etaval[neta] = eta;
                v1val[neta] = fabs(v1A);
                v1vale[neta] = fabs(v1Ae);
                neta++;
            }
            if (i == 0) gN1MCm22SUB2[cbin] = new TGraphErrors(netabins-i-1, deltaEta, v1val, 0, v1vale);
            if (i == 1) gN1MCm18SUB2[cbin] = new TGraphErrors(netabins-i-2, deltaEta, v1val, 0, v1vale);
            if (i == 2) gN1MCm14SUB2[cbin] = new TGraphErrors(netabins-i-3, deltaEta, v1val, 0, v1vale);
            if (i == 3) gN1MCm10SUB2[cbin] = new TGraphErrors(netabins-i-4, deltaEta, v1val, 0, v1vale);
            if (i == 4) gN1MCm06SUB2[cbin] = new TGraphErrors(netabins-i-5, deltaEta, v1val, 0, v1vale);
            if (i == 5) gN1MCm02SUB2[cbin] = new TGraphErrors(netabins-i-6, deltaEta, v1val, 0, v1vale);
            if (i == 6) gN1MCp02SUB2[cbin] = new TGraphErrors(netabins-i-6, deltaEta, v1val, 0, v1vale);
            if (i == 7) gN1MCp06SUB2[cbin] = new TGraphErrors(netabins-i-5, deltaEta, v1val, 0, v1vale);
            if (i == 8) gN1MCp10SUB2[cbin] = new TGraphErrors(netabins-i-4, deltaEta, v1val, 0, v1vale);
            if (i == 9) gN1MCp14SUB2[cbin] = new TGraphErrors(netabins-i-3, deltaEta, v1val, 0, v1vale);
            if (i == 10) gN1MCp18SUB2[cbin] = new TGraphErrors(netabins-i-2, deltaEta, v1val, 0, v1vale);
            if (i == 11) gN1MCp22SUB2[cbin] = new TGraphErrors(netabins-i-1, deltaEta, v1val, 0, v1vale);
            /*for (int j = 0; j<neta; j++) {
                cout<<"deltaEta = "<<deltaEta[j]<<"\tv1val = "<<v1val[j]<<"\tv1vale = "<<v1vale[j]<<endl;
            }*/
        }
    }

    // draw options
    for (int cbin = 0; cbin<cbinsCENT; cbin++) {
        gN1MCm22SUB2[cbin]->SetMarkerColor(kBlue);
        gN1MCm18SUB2[cbin]->SetMarkerColor(kRed);
        gN1MCm14SUB2[cbin]->SetMarkerColor(kGreen+2);
        gN1MCm10SUB2[cbin]->SetMarkerColor(kMagenta);
        gN1MCm06SUB2[cbin]->SetMarkerColor(kOrange+7);
        gN1MCm02SUB2[cbin]->SetMarkerColor(kBlack);
        gN1MCp02SUB2[cbin]->SetMarkerColor(kBlack);
        gN1MCp06SUB2[cbin]->SetMarkerColor(kOrange+7);
        gN1MCp10SUB2[cbin]->SetMarkerColor(kMagenta);
        gN1MCp14SUB2[cbin]->SetMarkerColor(kGreen+2);
        gN1MCp18SUB2[cbin]->SetMarkerColor(kRed);
        gN1MCp22SUB2[cbin]->SetMarkerColor(kBlue);

        gN1MCm22SUB2[cbin]->SetLineColor(kBlue);
        gN1MCm18SUB2[cbin]->SetLineColor(kRed);
        gN1MCm14SUB2[cbin]->SetLineColor(kGreen+2);
        gN1MCm10SUB2[cbin]->SetLineColor(kMagenta);
        gN1MCm06SUB2[cbin]->SetLineColor(kOrange+7);
        gN1MCm02SUB2[cbin]->SetLineColor(kBlack);
        gN1MCp02SUB2[cbin]->SetLineColor(kBlack);
        gN1MCp06SUB2[cbin]->SetLineColor(kOrange+7);
        gN1MCp10SUB2[cbin]->SetLineColor(kMagenta);
        gN1MCp14SUB2[cbin]->SetLineColor(kGreen+2);
        gN1MCp18SUB2[cbin]->SetLineColor(kRed);
        gN1MCp22SUB2[cbin]->SetLineColor(kBlue);

        gN1MCm22SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm18SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm14SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm10SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm06SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm02SUB2[cbin]->SetMarkerStyle(21);
        gN1MCp02SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp06SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp10SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp14SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp18SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp22SUB2[cbin]->SetMarkerStyle(25);

        gN1MCm22SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm18SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm14SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm10SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm06SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm02SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp02SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp06SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp10SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp14SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp18SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp22SUB2[cbin]->SetMarkerSize(1.2);
    }

    int setcent = 13;

    TCanvas * c0 = new TCanvas("c0","c0",650,600);
    TPad * pad0 = (TPad *) c0->cd();
    pad0->SetGrid();
    TH1D * h0 = new TH1D("h0", "h0", 100, 0, 4.4);
    h0->SetTitle("");
    h0->SetStats(0);
    h0->GetYaxis()->SetRangeUser(0, 0.06);
    h0->SetXTitle("#Delta#eta");
    h0->SetYTitle("|v_{1}^{even}|");
    h0->Draw();
    gN1MCm22SUB2[setcent-1]->Draw("same p");
    gN1MCm18SUB2[setcent-1]->Draw("same p");
    gN1MCm14SUB2[setcent-1]->Draw("same p");
    gN1MCm10SUB2[setcent-1]->Draw("same p");
    gN1MCm06SUB2[setcent-1]->Draw("same p");
    gN1MCm02SUB2[setcent-1]->Draw("same p");
    TArrow * ar1 = new TArrow(0.2, 0.2, 0.1, 0.7);
    ar1->SetAngle(40);
    ar1->SetLineWidth(2);
    ar1->SetLineColor(kBlue);
    ar1->Draw("same");
    TLegend * leg0 = new TLegend(0.61, 0.57, 0.88, 0.88);
    SetLegend(leg0,18);
    leg0->SetHeader(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    leg0->AddEntry(gN1MCm22SUB2[setcent-1],"N1MCm22SUB2","p");
    leg0->AddEntry(gN1MCm18SUB2[setcent-1],"N1MCm18SUB2","p");
    leg0->AddEntry(gN1MCm14SUB2[setcent-1],"N1MCm14SUB2","p");
    leg0->AddEntry(gN1MCm10SUB2[setcent-1],"N1MCm10SUB2","p");
    leg0->AddEntry(gN1MCm06SUB2[setcent-1],"N1MCm06SUB2","p");
    leg0->AddEntry(gN1MCm02SUB2[setcent-1],"N1MCm02SUB2","p");
    leg0->Draw();

}
