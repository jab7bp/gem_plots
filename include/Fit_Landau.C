void Fit_Landau(TH1F *histo){
  
  histo->Fit("landau","q");
  
  gPad->Update();
  TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
  st->SetOptStat(10);
      
  st->SetX1NDC(0.5);
  st->SetX2NDC(0.98);
  st->SetY1NDC(0.5);
  st->SetY2NDC(0.93);
  
}