// clusters on one layer
struct LayerCluster
{
  int layerID;
  ChamberCluster chamber_cluster[4]; // 4 chambers per layer
  TString GEMtype;

    LayerCluster():
        layerID(-1)
    {
    }

    void addCluster(const GEMCluster &c) 
    {
        int id = c.moduleID;
        if(id < 0 || id > 3) {
            cout<<"incorrect chamber pos index: "<< id<<" on layer: "<<layerID<<endl;
            return;
        }
        chamber_cluster[id].addCluster(c);
	chamber_cluster[id].moduleID = id;
    
    }

    void Match() {
      for(auto &i: chamber_cluster)
            i.Match();
    }

  void Clear(){
    for(auto &chamber: chamber_cluster)
     chamber.Hits2D.clear();
  }
};

