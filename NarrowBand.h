#ifndef NARROWBAND_H__
#define NARROWBAND_H__

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <utility>

#include<Eigen/Core>

#include <pcl/common/common.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkXMLImageDataWriter.h>

#include "colors.h"


typedef pcl::PointXYZ PointType;


class NarrowBand
{
public:
    NarrowBand(std::string filename);
    void build3Dgrid(int prec);
    void writeGridToFile();
    void computeInitialSurface();
    void writeUtoFile();
    void expandGrid(int count);
    void computeDistanceFunction();
    void writeDtoFile();
    void computeGradient();
    void writeDDtoFile();
    void evolve(int iter, float deltaT);
private:
    pcl::PointCloud<PointType>::Ptr _cloud;

    int _iter;
    double _res;
    float _delta;
    int _cDim,_vDim;
    int _numElem,_numCell;
    Eigen::Vector3f _min;

    std::vector<int> _cMask,_vMask;
    std::queue<int> _points;
    std::vector<float> _U,_Unew;
    std::vector<float> _D;
    std::vector<Eigen::Vector3f> _DD;

    float interp(float xp, float yp, float zp);
    double computeCloudResolution();
    float solve(std::vector<float> &a);
    int fromIJK2int(int i, int j, int k, int dim);
    void fromIJK2Center(int i, int j, int k, float out[]);
    void fromIJK2Vertex(int i, int j, int k, float out[]);
    void fromCenter2Cell(float x, float y, float z,int out[]);
    void fromint2IJK(int in, int dim, int out[]);
    void fromPoint2IJK(float x, float y, float z,int out[]);
    void fromVertex2IJK(float x, float y, float z,int out[]);

    int tagCell(int i, int j, int k);

};

class CompareDist
{
public:
    bool operator()(std::pair<int,float> n1,std::pair<int,float> n2) {
        return n1.second < n2.second;
    }
};

typedef std::priority_queue<std::pair<int,float>,std::vector<std::pair<int,float> >,CompareDist> PriorityQueue;


#endif //NARROWBAND_H__
