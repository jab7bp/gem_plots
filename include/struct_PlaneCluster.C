// clusters on one plane
struct PlaneCluster
{
    int axis;
    vector<GEMCluster> clusters;

    PlaneCluster():
        axis(-1)
    {
        clusters.clear();
    }

    void addCluster(const GEMCluster& c)
    {
        clusters.push_back(c);
    }

    // must be in desending order
    void Sort() {
        sort(clusters.begin(), clusters.end(), 
                [](const GEMCluster &c1, const GEMCluster &c2)
                {
                return c1.adc > c2.adc; // desending order
                });
    }

    void Clear(){
        clusters.clear();
    }
};

