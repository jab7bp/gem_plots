// show 2d cluster map for 1 layer
TH2F* GetLayerCluster2DMap(const LayerCluster &layer_cluster, const int &nlayer)
{
    TH2F *hClusterMap = new TH2F(Form("hClusterMap%d", nlayer), 
            Form("cluster 2d map layer %d", nlayer), 200, -400, 400, 600, -1200, 1200);
    for(auto &chamber : layer_cluster.chamber_cluster) {
        for(auto &c: chamber.Hits2D)
            hClusterMap -> Fill(c.xc.pos, c.yc.pos);
    }

    hClusterMap -> GetXaxis() -> SetTitle("X [mm]");
    hClusterMap -> GetXaxis() -> CenterTitle();
    hClusterMap -> GetYaxis() -> SetTitle("Y [mm]");
    hClusterMap -> GetYaxis() -> CenterTitle();
    hClusterMap -> GetYaxis() -> SetTitleOffset(0.8);

    return hClusterMap;
}