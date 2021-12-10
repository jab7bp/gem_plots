void Reset_Stats(TH1F *histo){
  
  gPad->Update();
  TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
      
  st->SetX1NDC(0.78);
  st->SetX2NDC(0.98);
  st->SetY1NDC(0.78);
  st->SetY2NDC(0.93);
  
}