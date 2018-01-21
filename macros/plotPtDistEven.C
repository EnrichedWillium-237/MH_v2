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
TString etags[] = {"-2.4", "-2.0", "-1.6", "-1.2", "-0.8", "-0.4", "00.0", "00.4", "00.8", "01.2", "01.6", "02.0", "02.4"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCm10SUB2", "N1MCm06SUB2", "N1MCm02SUB2",
                 "N1MCp02SUB2", "N1MCp06SUB2", "N1MCp10SUB2", "N1MCp14SUB2", "N1MCp18SUB2", "N1MCp22SUB2",
                 "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCm10SUB3", "N1MCm06SUB3", "N1MCm02SUB3",
                 "N1MCp02SUB3", "N1MCp06SUB3", "N1MCp10SUB3", "N1MCp14SUB3", "N1MCp18SUB3", "N1MCp22SUB3"
};

TH1D * gN1MCm22SUB3_pt[netabins][ncentbins];
TH1D * gN1MCm18SUB3_pt[netabins][ncentbins];
TH1D * gN1MCm14SUB3_pt[netabins][ncentbins];

TH1D * gN1MCp22SUB3_pt[netabins][ncentbins];
TH1D * gN1MCp18SUB3_pt[netabins][ncentbins];
TH1D * gN1MCp14SUB3_pt[netabins][ncentbins];

TH1D * gN1MC22SUB3_pt[netabins][ncentbins];
TH1D * gN1MC18SUB3_pt[netabins][ncentbins];
TH1D * gN1MC14SUB3_pt[netabins][ncentbins];

void plotPtDistEven()
{

    TH1::SetDefaultSumw2();

    for (int ebin = 0; ebin<netabins; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            gN1MCm22SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCm22SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCm18SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCm18SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCm14SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCm14SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);

            gN1MCp22SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCp22SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCp18SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCp18SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCp14SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCp14SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);

            gN1MC22SUB3_pt[ebin][cbin] = new TH1D(Form("N1MC22SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MC18SUB3_pt[ebin][cbin] = new TH1D(Form("N1MC18SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MC14SUB3_pt[ebin][cbin] = new TH1D(Form("N1MC14SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
        }
    }

    for (int i = 0; i<NANALS; i++) {
        for (int ebin = 0; ebin<netabins; ebin++) {
            for (int cbin = 0; cbin<ncentbins; cbin++) {
                TString tag = Form("figures_MH/%s/eta_%s_%s/data/%s_%d_%d_eta_%1.1f_%1.1f_A.dat",ANAL[i].Data(),etags[ebin].Data(),etags[ebin+1].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin],etabins[ebin],etabins[ebin+1]);
                cout<<tag<<endl;
            }
        }
    }
    // for (int cbin = 0; cbin<ncentbins; cbin++) {
    //     TString tag = Form("figures_MH/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
    //     // cout<<tag.Data()<<endl;
    //     ifstream fin(tag.Data());
    //
    //     double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
    //     int neta = 0;
    //     while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
    //     // if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;
    //
    //     if (i == 0) {
    //         gN1MCm22SUB2_eta[cbin]->SetBinContent(neta+1, v1A);
    //         gN1MCm22SUB2_eta[cbin]->SetBinError(neta+1, v1Ae);
    //     }
    // }

}
