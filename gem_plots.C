// Script to make various plots for GEM Modules and Layers
// Written/Modified/Compiled by John Boyd
// Last update: Dec. 24, 2021
// Update NOTES: fixed locatin of plots folder. Now make directories for pdf and png

// This script makes all the usual diagnositc plots for a set of GEM modules.
// To run simply check that gem_config.cfg matches your GEM setup. Also make sure that 
// the hit and cluster root files for this run are located in the ../Rootfiles directory
// and are named "cluster_bbgem_#runnumber.root" and "hit_bbgem_#runnumber.root".
//
// To run execute in terminal: root 'module_plots.C(#run number)'
//
//
// The following plots will be written to a pdf named ../plots/clusters/cluster_#runnumber.pdf
//
// X/Y strips hit
// X/Y strips ADC
// X/Y cluster size
// X/Y clusters per event
// X/Y clusters hit maps
// X/Y clusters ADC
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <sstream>
#include "./include/cluster_matching.C"
#include "./include/convert_uv_to_xy.C"
#include "./include/Rescale_histo.C"
#include "./include/Fit_Landau.C"
#include "./include/Reset_Stats.C"
#include "./include/set_styles.C"


#include "./include/ConfigParser.C"
#include "./include/struct_GEMCluster.C"
#include "./include/struct_GEM2DHit.C"
#include "./include/struct_PlaneCluster.C"
#include "./include/struct_ChamberCluster.C"
#include "./include/struct_LayerCluster.C"
#include "./include/TH2F_GetLayerCluster2DMap.C"

struct GEMCluster;
struct GEM2DHit;
struct PlaneCluster;
struct ChamberCluster;
struct LayerCluster;


void gem_plots(int run = 1000, bool plot_all = false, bool save = true, bool plot_cluster_maps = true, bool plot_strip_info = false, bool plot_clusters = false, bool plot_adc = false, bool nclust_size = false, int nmax = 100000000, bool skip_events = false){

  if(plot_all){
    plot_cluster_maps = true;
    plot_strip_info = true;
    plot_clusters = true;
    plot_adc = true;
    nclust_size = true;
  }

  char *user = getenv("USER");

  TString HitDir = "../../Rootfiles/";

  TString output = Form("../../plots/clusters/all_plots_%i.pdf",run);
  TString clust_output = Form("cluster_%i.pdf",run);

  gStyle->SetPalette(1);

  TChain *t_clust = new TChain("GEMCluster");
  TChain *t_hit = new TChain("GEMHit");

  t_clust->Add(HitDir + Form("cluster_0_*_%i.root",run));  //Cluster root file
  t_hit->Add(HitDir + Form("hit_*_%i.root",run));        //Hit root file

  ////Parse file name to use info for run labels on graphs ///
  TString t_clust_name = t_clust->GetFile()->GetName();
  TObjArray *token_tcn = t_clust_name.Tokenize("_");
  TObjArray *token_file_name = clust_output.Tokenize(".");
  TString prefix = token_tcn->At(2)->GetName();
  TString runInfo = prefix + " " + Form("%i", run);
  TString file_name = token_file_name->At(0)->GetName();
  
  TString clust_output_dir = "../../plots/clusters/" + prefix + "_" + run;
  mkdir(clust_output_dir,0777);

  ConfigParser("./include/gem_config.cfg");       //File that lists GEM layers and positions

  int nlayermax = 9;         //Should not have more layers or modules than this
  int nmodmax = 4;

  //Use if some events are prescaled
  bool skip_event = false;
  int remainder = 82;
  int modulo = 101;
  
  LayerCluster layer_cluster[nlayermax];
  int maxcluster = 5000;     //Events would never have more than 5000 clusters
  int maxch = 60000;   // Events never have more than 10000 channels

  int evtID;
  int nCluster;
  int layer[maxch];    
  int prodID[maxcluster];    
  int moduleID[maxch];
  int axis[maxch];
  int size[maxcluster];
  float adc[maxcluster];
  float pos[maxcluster];



  t_clust->SetBranchAddress("evtID",&evtID);
  t_clust->SetBranchAddress("nCluster",&nCluster);
  t_clust->SetBranchAddress("size",&size);
  t_clust->SetBranchAddress("planeID",layer); 
  t_clust->SetBranchAddress("axis",axis); 
  t_clust->SetBranchAddress("moduleID",moduleID); 
  t_clust->SetBranchAddress("adc",adc); 
  t_clust->SetBranchAddress("pos",pos); 
  

  //We define our histograms
  TH2F *cluster2D[nlayermax][nmodmax];
  TH2F *cluster2D_all[nlayermax];

  TH1F *hit_x[nlayermax][nmodmax];
  TH1F *hit_y[nlayermax][nmodmax];
  
  TH2F *ADC2D[nlayermax][nmodmax];
  TH1F *ADC_x[nlayermax][nmodmax];
  TH1F *ADC_y[nlayermax][nmodmax];

  TH1F *clust_size_x[nlayermax][nmodmax];
  TH1F *clust_size_y[nlayermax][nmodmax];
  TH1F *nclust_x[nlayermax][nmodmax];
  TH1F *nclust_y[nlayermax][nmodmax];
  
  TH1F *hit_x_strip[nlayermax][nmodmax];
  TH1F *hit_y_strip[nlayermax][nmodmax];
  TH1F *ADC_x_strip[nlayermax][nmodmax];
  TH1F *ADC_y_strip[nlayermax][nmodmax];

  double hit_strip_max[3] = {128*12, 128*30, 128*10};
  double clust_pos_L[nlayermax][nmodmax][2]; //0 is x and 1 is y
  double clust_pos_R[nlayermax][nmodmax][2];

  //Canvas for Strip info             
  TCanvas *c; //Do we need this canvas?
  if(plot_all) c = new TCanvas("c","",1600,1200);   
 
  //Canvas for Cluster size info
  TCanvas *c2;
  if(nclust_size) c2 = new TCanvas("c2","",1600,1200);

  //Canvas for Cluster pos info
  TCanvas *c3;
  if(plot_clusters) c3 = new TCanvas("c3","",1600,1200); 

  //Canvas for cluster ADC info
  TCanvas *c4;
  if(plot_adc) c4  = new TCanvas("c4","",1600,1200);

  //We initialize all the histograms
  int first_layer = 1;
  for(int ilayer=first_layer; ilayer<layer_list.size(); ilayer++){
  
    layer_cluster[layer_list[ilayer]].GEMtype = GEMtype[layer_list[ilayer]];
    layer_cluster[ilayer].layerID = ilayer;

    double side_offsets = 30;
    double left_x_all = -mod_Lx[layer_list[ilayer]]*(nmodules[layer_list[ilayer]]*1.0/2) - side_offsets;
    double right_x_all = -mod_Lx[layer_list[ilayer]]*(nmodules[layer_list[ilayer]]*1.0/2 - nmodules[layer_list[ilayer]]) + side_offsets;
    double left_y_all = -mod_Ly[layer_list[ilayer]]/2 - side_offsets;
    double right_y_all = mod_Ly[layer_list[ilayer]]/2 + side_offsets;

    int nbinsx = 800;
    

    if(GEMtype[layer_list[ilayer]] == "UV") nbinsx = 600;
	
    if(plot_cluster_maps) { 	      
            cluster2D_all[layer_list[ilayer]] = new TH2F(Form("cluster2D_all_%i",layer_list[ilayer]),Form("Cluster Distribution Map, Layer %i;x (mm);y (mm)",layer_list[ilayer]),nbinsx,left_x_all,right_x_all,200,left_y_all,right_y_all);
        }

    for(int imod = 0; imod < nmodules[layer_list[ilayer]]; imod++){

      
      double left_x = -mod_Lx[layer_list[ilayer]]*(nmodules[layer_list[ilayer]]*1.0/2 - imod) - side_offsets;
      double right_x = -mod_Lx[layer_list[ilayer]]*(nmodules[layer_list[ilayer]]*1.0/2 - imod - 1) + side_offsets;
      
      double left_y = -mod_Ly[layer_list[ilayer]]/2 - side_offsets;
      double right_y = mod_Ly[layer_list[ilayer]]/2 + side_offsets;

      if(GEMtype[layer_list[ilayer]] == "INFN"){
	double temp = left_x;
	left_x = -1*right_x;
	right_x = -1*temp;
      }

      clust_pos_L[layer_list[ilayer]][imod][0] = left_x;
      clust_pos_L[layer_list[ilayer]][imod][1] = left_y;
      clust_pos_R[layer_list[ilayer]][imod][0] = right_x;
      clust_pos_R[layer_list[ilayer]][imod][1] = right_y;
    
      
      double clust_y_left = left_y;
      double clust_y_right = right_y;

      if(GEMtype[layer_list[ilayer]] == "UV"){
	clust_y_left = left_x;
	clust_y_right = right_x;
      }

      int ADC_right = 2000;

      if(GEMtype[layer_list[ilayer]] == "UVa")
	ADC_right = 600;


      hit_x[layer_list[ilayer]][imod] = new TH1F(Form("hit_x_%i_%i",layer_list[ilayer],imod),Form("X Cluster Pos Layer %i Module %i;x (mm);",layer_list[ilayer],imod),200,left_x,right_x);
      hit_y[layer_list[ilayer]][imod] = new TH1F(Form("hit_y_%i_%i",layer_list[ilayer],imod),Form("Y Cluster Pos Layer %i Module %i;y (mm);",layer_list[ilayer],imod),200,clust_y_left,clust_y_right);
      cluster2D[layer_list[ilayer]][imod] = new TH2F(Form("cluster2D_%i_%i",layer_list[ilayer],imod),Form("Cluster Distribution Layer %i Module %i;x (mm);y (mm)",layer_list[ilayer],imod),nbinsx,left_x,right_x,200,left_y,right_y);

      ADC_x[layer_list[ilayer]][imod] = new TH1F(Form("ADC_x_%i_%i",layer_list[ilayer],imod),Form("X Cluster Charge Layer %i Module %i;ADC;",layer_list[ilayer],imod),200,0,ADC_right);
      ADC_y[layer_list[ilayer]][imod] = new TH1F(Form("ADC_y_%i_%i",layer_list[ilayer],imod),Form("Y Cluster Charge Layer %i Module %i;ADC;",layer_list[ilayer],imod),200,0,ADC_right);
      ADC2D[layer_list[ilayer]][imod] = new TH2F(Form("ADC2D_%i_%i",layer_list[ilayer],imod),Form("Cluster Charge Sharing Layer %i Module %i;ADC x;ADC y",layer_list[ilayer],imod),200,0,ADC_right,200,0,ADC_right);
      
      clust_size_x[layer_list[ilayer]][imod] = new TH1F(Form("clust_size_x_%i_%i",layer_list[ilayer],imod),Form("X Cluster Size  Layer %i Module %i;Size;",layer_list[ilayer],imod),15,0,15);
      clust_size_y[layer_list[ilayer]][imod] = new TH1F(Form("clust_size_y_%i_%i",layer_list[ilayer],imod),Form("Y Cluster Size Layer %i Module %i;Size;",layer_list[ilayer],imod),15,0,15);
      
      
      nclust_x[layer_list[ilayer]][imod] = new TH1F(Form("nclust_x_%i_%i",layer_list[ilayer],imod),Form("X Clusters per Event Layer %i Module %i;# of Clusters;",layer_list[ilayer],imod),30,0,30);
      nclust_y[layer_list[ilayer]][imod] = new TH1F(Form("nclust_y_%i_%i",layer_list[ilayer],imod),Form("Y Clusters per Event Layer %i Module %i;# of Clusters;",layer_list[ilayer],imod),30,0,30);
      
      if(plot_strip_info){
        ADC_x_strip[layer_list[ilayer]][imod] = new TH1F(Form("ADC_x_strip_%i_%i",layer_list[ilayer],imod),Form("X Strip ADC Layer %i Module %i;ADC;",layer_list[ilayer],imod),200,0,2000);
        ADC_y_strip[layer_list[ilayer]][imod] = new TH1F(Form("ADC_y_strip_%i_%i",layer_list[ilayer],imod),Form("Y Strip ADC Layer %i Module %i;ADC;",layer_list[ilayer],imod),200,0,2000);
        }   
      
      int nstrips = 0;

      if(GEMtype[layer_list[ilayer]] == "UVa") nstrips = hit_strip_max[0];
      if(GEMtype[layer_list[ilayer]] == "UV") nstrips = hit_strip_max[1];
      if(GEMtype[layer_list[ilayer]] == "INFN") nstrips = hit_strip_max[2];

      if(plot_strip_info){
        hit_x_strip[layer_list[ilayer]][imod] = new TH1F(Form("strip_hit_x_%i_%i",layer_list[ilayer],imod),Form("X Strip Pos Layer %i Module %i;strip;",layer_list[ilayer],imod),nstrips,0,nstrips);
        hit_y_strip[layer_list[ilayer]][imod] = new TH1F(Form("strip_hit_y_%i_%i",layer_list[ilayer],imod),Form("Y Strip Pos Layer %i Module %i;strip;",layer_list[ilayer],imod),nstrips, 0,nstrips);
        }

            
    }
  }
  
 
  vector<float> clust_x[nlayermax][nmodmax], clust_y[nlayermax][nmodmax];      //Fill variables with cluster positions
  vector<int> clust_x_i[nlayermax][nmodmax], clust_y_i[nlayermax][nmodmax];    //Fill variables with cluster index


  int ievent = 0;
  int n_ent_clust = t_clust->GetEntries();
  int n_ent_hit = t_hit->GetEntries();

  cout<<"Begin Clusters"<<endl;

  /////////////////// Get cluster information and fill histograms ////////////////////////////  
  while(t_clust->GetEntry(ievent++)){
    //cout<<ievent<<endl;
    if(ievent%10000==0) cout<<std::setprecision(2)<<ievent*1.0*100/n_ent_clust<<" % \r"<<std::flush;
    
    if(evtID > nmax) break;

    if(skip_event)
      if(evtID%modulo == remainder)
	continue;
     
    if(nCluster > 500) continue;   //Skip more than 50 clusters
    
    // Loop over all clusters in event
    for(int iclust = 0; iclust < nCluster; iclust++){
      
      GEMCluster cluster;
      cluster.layerID = layer[iclust];
      cluster.moduleID = moduleID[iclust];
      cluster.axis = axis[iclust];
      cluster.size = size[iclust];
      cluster.adc = adc[iclust];
      cluster.pos = pos[iclust];
     
      layer_cluster[layer[iclust]].addCluster(cluster);

    }

    // match clusters
    for(int i=0;i<nlayermax;i++){
         
      for(auto &chamber : layer_cluster[i].chamber_cluster){
	for(auto &plane : chamber.plane_cluster){

	  nCluster = plane.clusters.size();
	  if(nCluster > 0){
	    if(plane.axis == 0) nclust_x[layer_cluster[i].layerID][chamber.moduleID]->Fill(nCluster);
	    if(plane.axis == 1) nclust_y[layer_cluster[i].layerID][chamber.moduleID]->Fill(nCluster);
	  }

	  for(auto &cluster : plane.clusters){

	    if(cluster.axis == 0){
	      hit_x[cluster.layerID][cluster.moduleID]->Fill(cluster.pos);
	      ADC_x[cluster.layerID][cluster.moduleID]->Fill(cluster.adc);
	      clust_size_x[cluster.layerID][cluster.moduleID]->Fill(cluster.size);
	    }
	    if(cluster.axis == 1){
	      hit_y[cluster.layerID][cluster.moduleID]->Fill(cluster.pos);
	      ADC_y[cluster.layerID][cluster.moduleID]->Fill(cluster.adc);
	      clust_size_y[cluster.layerID][cluster.moduleID]->Fill(cluster.size);
	    }
	  }
	}
      }
      
      layer_cluster[i].Match();
      
      
      for(auto &chamber : layer_cluster[i].chamber_cluster){
	for(auto &hit2d: chamber.Hits2D){

	  if(layer_cluster[i].GEMtype == "UV")
	    convert_uv_to_xy(hit2d.xc.pos,hit2d.yc.pos,60);

      if(plot_clusters){
            cluster2D[hit2d.xc.layerID][hit2d.xc.moduleID]->Fill(hit2d.xc.pos,hit2d.yc.pos);
	       }

      if(plot_adc){
        ADC2D[hit2d.xc.layerID][hit2d.xc.moduleID]->Fill(hit2d.xc.adc,hit2d.yc.adc);
        }

      if(plot_cluster_maps)
        {
	       cluster2D_all[hit2d.xc.layerID]->Fill(hit2d.xc.pos,hit2d.yc.pos);
        }
	}
      }
       
     layer_cluster[i].Clear(); 
    }


  }

  
  
  cout<<"Finished Clustering"<<endl;
  cout<<"Now do strip hits"<<endl;

  int nch;
  int strip[maxch];
  int adc_all[6][maxch];

  t_hit->SetBranchAddress("evtID",&evtID);
  t_hit->SetBranchAddress("nch",&nch);
  t_hit->SetBranchAddress("axis",axis); 
  t_hit->SetBranchAddress("strip",strip); 
  t_hit->SetBranchAddress("planeID",layer); 
  t_hit->SetBranchAddress("moduleID",moduleID); 
  t_hit->SetBranchAddress("adc0",adc_all[0]); 
  t_hit->SetBranchAddress("adc1",adc_all[1]); 
  t_hit->SetBranchAddress("adc2",adc_all[2]); 
  t_hit->SetBranchAddress("adc3",adc_all[3]); 
  t_hit->SetBranchAddress("adc4",adc_all[4]); 
  t_hit->SetBranchAddress("adc5",adc_all[5]);  



  /////////////////////// Get strip information and fill histograms ////////////////////////////

  ievent = 0;

  
  while(t_hit->GetEntry(ievent++)){
    
    if(ievent%10000==0) cout<<std::setprecision(2)<<ievent*1.0*100/n_ent_hit<<" % \r"<<std::flush;    

    if(evtID > nmax) break;
    
    if(skip_event) 
      if(evtID%modulo == remainder) continue;

    if(nch > 600) continue;  //Skip events with more than 500 strips
    
    for(int istrip = 0; istrip < nch; istrip++){
   
      double adc_max = 0;

      //if(layer[istrip] == 1 && moduleID[istrip] == 1) 
    

      //Find max adc value
      for(int iadc = 0; iadc < 6; iadc++)
	if(adc_all[iadc][istrip] > adc_max) adc_max = adc_all[iadc][istrip];
      
      //Fill histograms with ADC
    if(plot_strip_info)
        {
          if(axis[istrip] == 0){
    	       ADC_x_strip[layer[istrip]][moduleID[istrip]]->Fill(adc_max);
    	       hit_x_strip[layer[istrip]][moduleID[istrip]]->Fill(strip[istrip]);
            }
          if(axis[istrip] == 1){
    	       ADC_y_strip[layer[istrip]][moduleID[istrip]]->Fill(adc_max);
    	       hit_y_strip[layer[istrip]][moduleID[istrip]]->Fill(strip[istrip]);
    	//if(layer[istrip] == 4 && moduleID[istrip] == 1 && strip[istrip] > 60 && strip[istrip] < 100 ) cout<<evtID<<endl;
            }
        }
    }

  }



  //Canvas for all 2D hitmaps
  int ncan = layer_list.size()/3 + 1;
  TCanvas *c5[ncan];
  if(plot_cluster_maps)
    {
      for(int ican = 0; ican < ncan-1; ican++){
        c5[ican] = new TCanvas(Form("c5_%i",ican),"",1600,1200);
        c5[ican]->Divide(1,3);
        }
    }

  gStyle->SetOptFit(011);

  
  int count = 1;
  int first_canvas = 0;
  for(int ilayer = first_layer; ilayer<layer_list.size();ilayer++){
    
    int ican = ilayer/3;
    int c_i = ilayer - ican*3;
    int canvas_i = first_canvas/3;

    if(plot_cluster_maps)
        {
            c5[canvas_i]->cd(count);
            Rescale_histo(cluster2D_all[layer_list[ilayer]],clust_pos_L[layer_list[ilayer]][0][1],clust_pos_R[layer_list[ilayer]][nmodules[layer_list[ilayer]]-1][1]);
            cluster2D_all[layer_list[ilayer]]->Draw("colz");
            TPaveText *runLabel = new TPaveText(.87, .95, .97, .999, "NDC");
            runLabel->SetBorderSize(1);
            runLabel->AddText("Run info: " + runInfo);
            runLabel->Draw();
            gPad->Update();
        }

    
    //LABELS and LINES and ETC.
    if(plot_cluster_maps)
        {
           
            TPaveStats *st = (TPaveStats*)cluster2D_all[layer_list[ilayer]]->FindObject("stats");
            st->SetOptStat(11);
            
            st->SetX1NDC(0.85);
            st->SetX2NDC(0.98);
            st->SetY1NDC(0.85);
            st->SetY2NDC(0.93);

            double gem_label_0_x1 = .1875;
            double gem_label_0_x2 = .22;
            double gem_label_y1 = 0.06;
            double gem_label_y2 = 0.1;

            TPaveText *gemLabel0 = new TPaveText(gem_label_0_x1, gem_label_y1, gem_label_0_x2, gem_label_y2, "NDC");
            TPaveText *gemLabel1 = new TPaveText(gem_label_0_x1 + .2, gem_label_y1, gem_label_0_x2 + .2, gem_label_y2, "NDC");
            TPaveText *gemLabel2 = new TPaveText(gem_label_0_x1 + .4, gem_label_y1, gem_label_0_x2 + .4, gem_label_y2, "NDC");
            TPaveText *gemLabel3 = new TPaveText(gem_label_0_x1 + .6, gem_label_y1, gem_label_0_x2 + .6, gem_label_y2, "NDC");
            gemLabel0->AddText("GEM 0");
            gemLabel0->Draw();
            gemLabel1->AddText("GEM 1");
            gemLabel1->Draw();
            gemLabel2->AddText("GEM 2");
            gemLabel2->Draw();
            gemLabel3->AddText("GEM 3");
            gemLabel3->Draw();

            bool drawGemLines = true;
            if(drawGemLines)
            {
            TLine *gem00 = new TLine(-2*uva_xy_x, -uva_xy_y/2 - 40, -2*uva_xy_x, uva_xy_y/2 + 40);
            TLine *gem01 = new TLine(-uva_xy_x, -uva_xy_y/2 - 40, -uva_xy_x, uva_xy_y/2 + 40);
            TLine *gem12 = new TLine(0, -uva_xy_y/2 - 40, 0, uva_xy_y/2 + 40);
            TLine *gem23 = new TLine(uva_xy_x, -uva_xy_y/2 - 40, uva_xy_x, uva_xy_y/2 + 40);
            TLine *gem33 = new TLine(2*uva_xy_x, -uva_xy_y/2 - 40, 2*uva_xy_x, uva_xy_y/2 + 40);

            gemLines(gem00);
            gemLines(gem01);
            gemLines(gem12);
            gemLines(gem23);
            gemLines(gem33);
            }
            


        }
      count++;
      if(count == 4)
        {
          count = 1;
        }
      first_canvas++;
  }
  

  // Loop over all layesr in config file
  for(int ilayer=first_layer; ilayer<layer_list.size(); ilayer++){  
    
    //c->SetWindowSize(400*nmodules[layer_list[ilayer]],300*nmodules[layer_list[ilayer]]);
    //c2->SetCanvasSize(400*nmodules[layer_list[ilayer]],1200);
    //c3->SetCanvasSize(400*nmodules[layer_list[ilayer]],1200);
    //c4->SetCanvasSize(400*nmodules[layer_list[ilayer]],1200);
    if(plot_all) c->Divide(nmodules[layer_list[ilayer]],4);    
    if(nclust_size) c2->Divide(nmodules[layer_list[ilayer]],4);
    if(plot_clusters) c3->Divide(nmodules[layer_list[ilayer]],3);
    if(plot_adc) c4->Divide(nmodules[layer_list[ilayer]],3);

    //Loop over all modules in layer
    for(int imod=0; imod<nmodules[layer_list[ilayer]]; imod++){

      int nstrips = 0;
      if(GEMtype[layer_list[ilayer]] == "UVa") nstrips = hit_strip_max[0];
      if(GEMtype[layer_list[ilayer]] == "UV") nstrips = hit_strip_max[1];
      if(GEMtype[layer_list[ilayer]] == "INFN") nstrips = hit_strip_max[2];

      if(plot_strip_info)
        {
          //Canvas 1 with strip info
          c->cd(imod + 1);
          Rescale_histo(hit_x_strip[layer_list[ilayer]][imod],0,nstrips);
          hit_x_strip[layer_list[ilayer]][imod]->Draw();
          Reset_Stats(hit_x_strip[layer_list[ilayer]][imod]);

          c->cd(imod + 1 + nmodules[layer_list[ilayer]]*1);
          Rescale_histo(hit_y_strip[layer_list[ilayer]][imod],0,nstrips);
          hit_y_strip[layer_list[ilayer]][imod]->Draw();
          Reset_Stats(hit_y_strip[layer_list[ilayer]][imod]);

          c->cd(imod + 1 + nmodules[layer_list[ilayer]]*2);
          ADC_x_strip[layer_list[ilayer]][imod]->Draw();
          Fit_Landau(ADC_x_strip[layer_list[ilayer]][imod]);
          
          c->cd(imod + 1 + nmodules[layer_list[ilayer]]*3);
          ADC_y_strip[layer_list[ilayer]][imod]->Draw();
          Fit_Landau(ADC_y_strip[layer_list[ilayer]][imod]);
        }


      //Canvas 2 with cluster size info
        if(nclust_size)
            {
              c2->cd(imod + 1);	
            

              clust_size_x[layer_list[ilayer]][imod]->Draw();
              Reset_Stats(clust_size_x[layer_list[ilayer]][imod]);

              c2->cd(imod + 1 + nmodules[layer_list[ilayer]]*1);
              clust_size_y[layer_list[ilayer]][imod]->Draw();
              Reset_Stats(clust_size_y[layer_list[ilayer]][imod]);

              c2->cd(imod + 1 + nmodules[layer_list[ilayer]]*2);
              nclust_x[layer_list[ilayer]][imod]->Draw();
              Reset_Stats(nclust_x[layer_list[ilayer]][imod]);

              c2->cd(imod + 1 + nmodules[layer_list[ilayer]]*3);
              nclust_y[layer_list[ilayer]][imod]->Draw();
              Reset_Stats(nclust_y[layer_list[ilayer]][imod]);
            }

      //Canvas 3 with cluster pos info
      if(plot_clusters)
        {
          c3->cd(imod + 1);
          Rescale_histo(cluster2D[layer_list[ilayer]][imod],clust_pos_L[layer_list[ilayer]][imod][1],clust_pos_R[layer_list[ilayer]][imod][1]);
          cluster2D[layer_list[ilayer]][imod]->Draw("colz");
      

          //Cosmetic adjustments of stats box
          gPad->Update();
          TPaveStats *st = (TPaveStats*)cluster2D[layer_list[ilayer]][imod]->FindObject("stats");
          st->SetOptStat(11);

          st->SetX1NDC(0.75);
          st->SetX2NDC(0.98);
          st->SetY1NDC(0.85);
          st->SetY2NDC(0.93);

          c3->cd(imod + 1 + nmodules[layer_list[ilayer]]*1);
          Rescale_histo(hit_x[layer_list[ilayer]][imod],clust_pos_L[layer_list[ilayer]][imod][0],clust_pos_R[layer_list[ilayer]][imod][0]);
          hit_x[layer_list[ilayer]][imod]->Draw();
          Reset_Stats(hit_x[layer_list[ilayer]][imod]);

          c3->cd(imod + 1 + nmodules[layer_list[ilayer]]*2);
          Rescale_histo(hit_y[layer_list[ilayer]][imod],clust_pos_L[layer_list[ilayer]][imod][1],clust_pos_R[layer_list[ilayer]][imod][1]);
          hit_y[layer_list[ilayer]][imod]->Draw();
          Reset_Stats(hit_y[layer_list[ilayer]][imod]);
        }


      //Canvas 4 with cluster ADC info
        if(plot_adc)
            {
              c4->cd(imod + 1);
              Rescale_histo(ADC2D[layer_list[ilayer]][imod]);
              ADC2D[layer_list[ilayer]][imod]->Draw("colz");

              TLine *line = new TLine(0,0,1600,1600);
              line->SetLineColor(kRed);
              line->SetLineWidth(2);
              line->Draw("same");


              c4->cd(imod + 1 + nmodules[layer_list[ilayer]]*1);
              ADC_x[layer_list[ilayer]][imod]->Draw();
              Fit_Landau(ADC_x[layer_list[ilayer]][imod]);

              c4->cd(imod + 1 + nmodules[layer_list[ilayer]]*2);
              ADC_y[layer_list[ilayer]][imod]->Draw();
              Fit_Landau(ADC_y[layer_list[ilayer]][imod]);
            }

    }

    //Print canvas to pdf for this layer
    
    if(nclust_size) {
        c2->Print(output);
        c2->Clear();
        }

    if(plot_clusters) {
        c3->Print(output);
        c3->Clear();
        }

    if(plot_adc) {
        c4->Print(output);
        c4->Clear();
        }
    
  }

  
  //PRINT Cluster Maps
  TString pdf_dir = "../../plots/clusters/" + prefix + "_" + run + "/pdf";
  TString png_dir = "../../plots/clusters/" + prefix + "_" + run + "/png";

  if(save) 
  {
    mkdir(pdf_dir, 0777);
    if(plot_cluster_maps) 
    {
      for(int ican = 0; ican < ncan-1; ican++) {
        if(ican == 0)
        {
          c5[0]->Print(pdf_dir + "/" + clust_output + "("); //Creating the pdf document
        }
        if(ican == ncan - 2)
        {
          c5[ican]->Print(pdf_dir + "/" + clust_output + ")"); //Closing the pdf after last canvas
        }
        else  if(ican > 0) c5[ican]->Print(pdf_dir + "/" + clust_output); //Middle canvases
      }
    }
  }

  if(save) 
  {
    mkdir(png_dir, 0777);
    for(int ican = 0; ican < ncan-1; ican++) {
      if(plot_cluster_maps) 
      {
        c5[ican]->SaveAs("../../plots/clusters/" + prefix + "_" + run + "/" + "png/" + file_name + Form("_%i.png", ican));
      }
    }
  }


}
