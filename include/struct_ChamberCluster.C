// clusters on one chamber
struct ChamberCluster
{
    int moduleID;
    PlaneCluster plane_cluster[2]; // 2 planes: x/y

    vector<GEM2DHit> Hits2D;

    ChamberCluster():
        moduleID(-1)
    {
    }

    void addCluster(const GEMCluster &c)
    {
        int id = c.axis;
        if(id < 0 || id > 1) {
            cout<<"axis index is incorrect: "<<id<<endl;
            return;
        }
        plane_cluster[id].addCluster(c);
	plane_cluster[id].axis = id;
    }

    void Sort() {
        plane_cluster[0].Sort();
        plane_cluster[1].Sort();
    }

    void Match() 
    {
        Sort();

        // get 2d hits
        int nCluster = plane_cluster[0].clusters.size();
        if(nCluster > plane_cluster[1].clusters.size())
            nCluster = plane_cluster[1].clusters.size();

        for(int i=0;i<nCluster;i++) {
            GEM2DHit hit(plane_cluster[0].clusters[i],
                    plane_cluster[1].clusters[i]);
            Hits2D.push_back(hit);
        }

        // clear this event, get ready for next event
        for(auto &i: plane_cluster)
	  i.Clear();
    }
  
};

