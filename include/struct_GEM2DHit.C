// reconstructed gem 2d hit point
struct GEM2DHit 
{
    GEMCluster xc, yc;
    GEM2DHit()
    {}

    GEM2DHit(const GEMCluster &a, const GEMCluster &b)
        :xc(a), yc(b)
    {}

};

