void Rescale_histo(TH1F *histo, double left, double right){

  int dx = (right - left)/3;

  double max = 1e10;

  for(int i=0; i<3; i++){
    int bin_L = histo->GetXaxis()->FindBin(left + dx*i);
    int bin_R = histo->GetXaxis()->FindBin(left + dx*(i+1));

    histo->GetXaxis()->SetRange(bin_L,bin_R);
    double mval = histo->GetMaximum();

    if(mval < max) max = mval;
  }
  
					  

  histo->GetXaxis()->SetRange(0,histo->GetXaxis()->FindBin(right));
  histo->GetYaxis()->SetRangeUser(0,max*1.5);

}

void Rescale_histo(TH2F *histo, double left, double right){

  int dx = (right - left)/3;

  double max = 1e10;

  for(int i=0; i<3; i++){
    int bin_L = histo->GetYaxis()->FindBin(left + dx*i);
    int bin_R = histo->GetYaxis()->FindBin(left + dx*(i+1));

    histo->GetYaxis()->SetRange(bin_L,bin_R);
    double mval = histo->GetMaximum();

    if(mval < max) max = mval;
  }
  

  histo->GetYaxis()->SetRange(0,histo->GetYaxis()->FindBin(right));
  histo->GetZaxis()->SetRangeUser(0,max);

}

void Rescale_histo(TH2F *histo){

  double mval = histo->GetMaximum();
 
  histo->GetZaxis()->SetRangeUser(0,mval/2);

}