void flip2D(TH2D * &h) {
  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();
  for(int j = 1; j<=ny/2; j++) {
    for(int i = 1; i<=nx; i++) {
      double hold = h->GetBinContent(i,j);
      double holde = h->GetBinError(i,j);
      h->SetBinContent(i,j,h->GetBinContent(i,ny+1-j));
      h->SetBinError(i,j,h->GetBinError(i,ny+1-j));
      h->SetBinContent(i,ny+1-j,hold);
      h->SetBinError(i,ny+1-j,holde);
    }
  } 
}
void test(){
  double etabins[]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4};
  double ptbins[]={0.3,0.4,0.5,0.6,0.8,1.0};
  TH2D * h = new TH2D("h","h",5,ptbins,12,etabins);
  h->SetBinContent(2,2,5);
  h->SetBinError(2,2,.2);
  h->SetOption("colz");
  TH2D * hc = (TH2D *) h->Clone("orignal");
  TCanvas * c = new TCanvas("c","c",1000,800);
  c->Divide(2);
  c->cd(1);
  hc->Draw();
  flip2D(h);
  c->cd(2);
  h->Draw();
}
