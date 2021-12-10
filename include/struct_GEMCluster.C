// gem cluster struct
//file: GEMCluster.C
struct GEMCluster
{
    int layerID, detID, moduleID, axis, size;
    float adc, pos;
 

    GEMCluster():
        layerID(-1), detID(-1), moduleID(-1), axis(-1), size(-1),
        adc(0), pos(-9999)
    {
    }

    GEMCluster(const GEMCluster &that):
        layerID(that.layerID), detID(that.detID), moduleID(that.moduleID),
        axis(that.axis), size(that.size), adc(that.adc), pos(that.pos)
    {
    }

    GEMCluster(GEMCluster &&that):
        layerID(that.layerID), detID(that.detID), moduleID(that.moduleID),
        axis(that.axis), size(that.size), adc(that.adc), pos(that.pos)
    {
    }

    // copy assignment
    GEMCluster& operator=(const GEMCluster &rhs)
    {
        if(this == &rhs)
            return *this;
        GEMCluster c(rhs);
        *this = move(c);
        return *this;
    }

    // move assignment
    GEMCluster & operator=(GEMCluster &&rhs) {
        if(this == &rhs)
            return *this;
        layerID = rhs.layerID;
        detID = rhs.detID;
        moduleID = rhs.moduleID;
        axis = rhs.axis;
        size = rhs.size;
        adc = rhs.adc;
        pos = rhs.pos;
        return *this;
    }
};

