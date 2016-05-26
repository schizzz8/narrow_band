#include "NarrowBand.h"

int getdigit(double number, int digit) {
    div_t divresult = div(number/pow(10, digit),10);
    return  divresult.rem;
}

void manageDirectories(std::string directory);
int isDirectoryEmpty(const char *dirname);


NarrowBand::NarrowBand(std::string filename):_cloud(new pcl::PointCloud<PointType>),_iter(0) {
    pcl::io::loadPCDFile<PointType>(filename.c_str(),*_cloud);
    std::cout << "\n";
    std::cout << BOLD(FBLU("Loading the Data-Set:\n"));
    std::cout << "\t>> Filename: " << filename << "\n";
    std::cout << "\t>> Points: " << _cloud->size() << "\n";

    size_t lastindex = filename.find_last_of(".");
    std::string directory = filename.substr(0,lastindex);
    manageDirectories(directory);
}

void NarrowBand::build3Dgrid(int prec) {
    PointType min_pt,max_pt;
    pcl::getMinMax3D(*_cloud,min_pt,max_pt);
    for(int i = 0; i < _cloud->size(); i++)
        _cloud->at(i).x -= min_pt.x,_cloud->at(i).y -= min_pt.y,_cloud->at(i).z -= min_pt.z;
    pcl::getMinMax3D(*_cloud,min_pt,max_pt);
    for(int i = 0; i < _cloud->size(); i++)
        _cloud->at(i).x -= max_pt.x/2,_cloud->at(i).y -= max_pt.y/2,_cloud->at(i).z -= max_pt.z/2;

    pcl::getMinMax3D(*_cloud,min_pt,max_pt);
    _res = computeCloudResolution();
    std::cout << "\t>> Min: (" << min_pt.x << "," << min_pt.y << "," << min_pt.z << ")\n";
    std::cout << "\t>> Max: (" << max_pt.x << "," << max_pt.y << "," << max_pt.z << ")\n";
    std::cout << "\t>> Average distance: " << _res << "m\n";
    std::cout << "--------------------------------------------------------------------------------\n";

    pcl::PCLPointCloud2 out;
    pcl::toPCLPointCloud2(*_cloud,out);
    pcl::io::saveVTKFile("data_set.vtk",out);

    bool found = false;
    int pos = 2;
    while(found == false) {
        pos--;
        if(getdigit(_res,pos) != 0)
            found = true;
    }

    _delta = (100.0/prec)*pow(10,pos);

    std::vector<float> ranges;
    ranges.push_back(max_pt.x-min_pt.x);
    ranges.push_back(max_pt.y-min_pt.y);
    ranges.push_back(max_pt.z-min_pt.z);

    float m = -1000000000;
    for(int i = 0; i < ranges.size(); i++)
        if(ranges[i] > m)
            m = ranges[i];

    _vDim = m/_delta;
    if(_vDim%2 == 0)
        _vDim += 21;
    else
        _vDim += 20;

    _cDim = _vDim - 1;

    _numElem = _vDim*_vDim*_vDim;
    _numCell = _cDim*_cDim*_cDim;

    _cMask = std::vector<int> (_numCell,0);
    _vMask = std::vector<int> (_numElem,0);
    _U = std::vector<float> (_numElem,1);
    _Unew = std::vector<float> (_numElem,0);
    _D = std::vector<float> (_numElem,10);
    _DD = std::vector<Eigen::Vector3f> (_numElem);

    _min = Eigen::Vector3f (-_cDim/2*_delta,-_cDim/2*_delta,-_cDim/2*_delta);

    for(int ii=0; ii < _cloud->size(); ii++) {
        float x = _cloud->at(ii).x, y = _cloud->at(ii).y, z = _cloud->at(ii).z;
        int index[3];
        fromPoint2IJK(x,y,z,index);
        int ijk = fromIJK2int(index[0],index[1],index[2],_cDim);
        if(ijk >= 0 && ijk <= _numCell-1)
            _cMask[ijk] = 1;
        else
            std::cerr << "Index exceeds grid dimensions!\n";

        for(int c = 0; c < 2; c ++)
            for(int b = 0; b < 2; b ++)
                for(int a = 0; a < 2; a ++) {
                    float vertex[3];
                    fromIJK2Vertex(index[0]+a,index[1]+b,index[2]+c,vertex);
                    int idx = fromIJK2int(index[0]+a,index[1]+b,index[2]+c,_vDim);
                    float dist = sqrt(pow(x - vertex[0],2) + pow(y - vertex[1],2) + pow(z - vertex[2],2));
                    _D[idx] = std::min(dist,_D[idx]);
                    _vMask[idx] = 1;
                }
    }

    for(int i = 0; i < _numCell; i++)
        if(_cMask[i] == 1)
            _points.push(i);
}

void NarrowBand::expandGrid(int count) {
    int last = _points.back();
    while(count > 0) {
        int cGrid[3];
        int current = _points.front();
        fromint2IJK(current,_cDim,cGrid);
        for(int c = -1; c < 2; c++)
            for(int b = -1; b < 2; b++)
                for(int a = -1; a < 2; a++)
                    if(a!=0 || b!=0 || c!=0) {
                        int ijk = fromIJK2int(cGrid[0]+a,cGrid[1]+b,cGrid[2]+c,_cDim);
                        if(_cMask[ijk] == 0) {
                            _cMask[ijk] = 1;
                            _points.push(ijk);
                            for(int cc = 0; cc < 2; cc++)
                                for(int bb = 0; bb < 2; bb++)
                                    for(int aa = 0; aa < 2; aa++) {
                                        int idx = fromIJK2int(cGrid[0]+a+aa,cGrid[1]+b+bb,cGrid[2]+c+cc,_vDim);
                                        if(_vMask[idx] == 0) {
                                            _vMask[idx] = 1;
                                            _D[idx] = 2;
                                        }
                                    }
                        }
                    }
        if(current == last) {
            count--;
            last = _points.back();
        }
        _points.pop();
    }
}

void NarrowBand::writeGridToFile() {
    std::cout << "\n";
    std::cout << BOLD(FBLU("Building the 3D Grid:\n"));
    std::cout << "\t>> Cells per side (total): " << _cDim << " (" << _numCell << ")\n";
    std::cout << "\t>> Delta: " << _delta << "m\n";
    std::cout << "\t>> Origin: (" << _min(0) << "," << _min(1) << "," << _min(2) << ")\n";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(int k = 0; k < _vDim; k++)
        for(int j = 0; j < _vDim; j++)
            for(int i = 0; i < _vDim; i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    float world[3];
                    fromIJK2Vertex(i,j,k,world);
                    points->InsertNextPoint (world[0],world[1],world[2]);
                }
            }

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("vGrid.vtp");

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polydata);
#else
    writer->SetInputData(polydata);
#endif
    writer->SetDataModeToAscii();
    writer->Write();

    int count = 0;
    vtkSmartPointer<vtkPoints> cpoints = vtkSmartPointer<vtkPoints>::New();
    for(int k = 0; k < _cDim; k++)
        for(int j = 0; j < _cDim; j++)
            for(int i = 0; i < _cDim; i++) {
                int ijk = fromIJK2int(i,j,k,_cDim);
                if(_cMask[ijk] == 1) {
                    count++;
                    float world[3];
                    fromIJK2Center(i,j,k,world);
                    cpoints->InsertNextPoint (world[0],world[1],world[2]);
                }
            }

    std::cout << "\t>> Active Cells: " << count << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";


    vtkSmartPointer<vtkPolyData> cpolydata = vtkSmartPointer<vtkPolyData>::New();
    cpolydata->SetPoints(cpoints);

    vtkSmartPointer<vtkXMLPolyDataWriter> cwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    cwriter->SetFileName("cGrid.vtp");

#if VTK_MAJOR_VERSION <= 5
    cwriter->SetInput(cpolydata);
#else
    cwriter->SetInputData(cpolydata);
#endif
    cwriter->SetDataModeToAscii();
    cwriter->Write();
}

void NarrowBand::computeDistanceFunction() {
    std::cout << "\n";
    std::cout << BOLD(FBLU("Level-Set Evolution:\n"));
    std::clock_t t0 = clock();
    //#####
    //# 1 #
    //#####
    for(int k = 0; k < _vDim; k++)
        for(int j = 0; j < _vDim; j++)
            for(int i = 0; i < _vDim; i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################

    //#####
    //# 2 #
    //#####
    for(int k = 0; k < _vDim; k++)
        for(int j = 0; j < _vDim; j++)
            for(int i = _vDim-1; i >= 0; i--) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################

    //#####
    //# 3 #
    //#####
    for(int k = 0; k < _vDim; k++)
        for(int j = _vDim-1; j >= 0; j--)
            for(int i = _vDim-1; i >= 0; i--) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################

    //#####
    //# 4 #
    //#####
    for(int k = 0; k < _vDim; k++)
        for(int j = _vDim-1; j >= 0; j--)
            for(int i = 0; i < _vDim; i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################

    //#####
    //# 5 #
    //#####
    for(int k = _vDim-1; k >= 0; k--)
        for(int j = 0; j < _vDim; j++)
            for(int i = 0; i < _vDim; i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################

    //#####
    //# 6 #
    //#####
    for(int k = _vDim-1; k >= 0; k--)
        for(int j = 0; j < _vDim; j++)
            for(int i = _vDim-1; i >= 0; i--) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################

    //#####
    //# 7 #
    //#####
    for(int k = _vDim-1; k >= 0; k--)
        for(int j = _vDim-1; j >= 0; j--)
            for(int i = _vDim-1; i >= 0; i--) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################

    //#####
    //# 8 #
    //#####
    for(int k = _vDim-1; k >= 0; k--)
        for(int j = _vDim-1; j >= 0; j--)
            for(int i = 0; i < _vDim; i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    std::vector<float> D_min;
                    D_min.push_back( std::min(_D[im1jk],_D[ip1jk]));
                    D_min.push_back( std::min(_D[ijm1k],_D[ijp1k]));
                    D_min.push_back( std::min(_D[ijkm1],_D[ijkp1]));
                    std::sort(D_min.begin(),D_min.end());
                    _D[ijk]=std::min(solve(D_min),_D[ijk]);

                }
            }
    //#########################################################################


    std::clock_t t1 = clock();
    double elapsed_time1 = double(t1 - t0)/CLOCKS_PER_SEC;
    std::cout << "\tTime to compute Distance Function: " << elapsed_time1 << "s\n";
}

void NarrowBand::writeDtoFile() {
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(_vDim,_vDim,_vDim);
    imageData->SetOrigin(_min.x(),_min.y(),_min.z());
    imageData->SetSpacing(_delta,_delta,_delta);
    imageData->SetExtent(0,_vDim-1,0,_vDim-1,0,_vDim-1);
#if VTK_MAJOR_VERSION <= 5
    imageData->SetNumberOfScalarComponents(1);
    imageData->SetScalarTypeToFloat();
#else
    imageData->AllocateScalars(VTK_FLOAT, 1);
#endif

    vtkSmartPointer<vtkFloatArray> dist = vtkSmartPointer<vtkFloatArray>::New();
    dist->SetName("distance_function");

    for(int k = 0; k < _vDim; k++)
        for(int j = 0; j < _vDim; j++)
            for(int i = 0; i < _vDim; i++)
                dist->InsertNextValue(_D[fromIJK2int(i,j,k,_vDim)]);

    imageData->GetPointData()->AddArray(dist);
    //imageData->GetCellData()->AddArray(dist);
    imageData->Update();


    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName("distance_function.vti");
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();
}

void NarrowBand::computeGradient() {
    std::clock_t t0 = clock();
    for(int k=0; k < _vDim; k++)
        for(int j=0; j < _vDim; j++)
            for(int i=0; i < _vDim; i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                if(_vMask[ijk] == 1) {
                    float d_x,d_y,d_z;
                    int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                    int im1jk = fromIJK2int(i-1,j,k,_vDim);
                    int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                    int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                    int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                    int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                    if(_vMask[ip1jk] == 1 && _vMask[im1jk] == 1)
                        d_x = (_D[ip1jk] - _D[im1jk])/(2*_delta);
                    else
                        d_x = std::min((_D[ip1jk]-_D[ijk])/_delta,(_D[ijk]-_D[im1jk])/_delta);

                    if(_vMask[ijp1k] == 1 && _vMask[ijm1k] == 1)
                        d_y = (_D[ijp1k] - _D[ijm1k])/(2*_delta);
                    else
                        d_y = std::min((_D[ijp1k]-_D[ijk])/_delta,(_D[ijk]-_D[ijm1k])/_delta);

                    if(_vMask[ijkp1] == 1 && _vMask[ijkm1] == 1)
                        d_z = (_D[ijkp1] - _D[ijkm1])/(2*_delta);
                    else
                        d_z = std::min((_D[ijkp1]-_D[ijk])/_delta,(_D[ijk]-_D[ijkm1])/_delta);

                    _DD[ijk] = Eigen::Vector3f(d_x,d_y,d_z);
                }
            }
    std::clock_t t1 = clock();
    double elapsed_time1 = double(t1 - t0)/CLOCKS_PER_SEC;
    std::cout << "\tTime to compute Distance Function Gradient: " << elapsed_time1 << "s\n";
}

void NarrowBand::writeDDtoFile() {
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(_vDim,_vDim,_vDim);
    imageData->SetOrigin(_min.x(),_min.y(),_min.z());
    imageData->SetSpacing(_delta,_delta,_delta);
    imageData->SetExtent(0,_vDim-1,0,_vDim-1,0,_vDim-1);
#if VTK_MAJOR_VERSION <= 5
    imageData->SetNumberOfScalarComponents(1);
    imageData->SetScalarTypeToFloat();
#else
    imageData->AllocateScalars(VTK_FLOAT, 1);
#endif

    vtkSmartPointer<vtkFloatArray> grad = vtkSmartPointer<vtkFloatArray>::New();
    grad->SetName("gradient");
    grad->SetNumberOfComponents(3);

    for(int k = 0; k < _vDim; k++)
        for(int j = 0; j < _vDim; j++)
            for(int i = 0; i < _vDim; i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                float tuple[3] = {_DD[ijk].x(),_DD[ijk].y(),_DD[ijk].z()};
                grad->InsertNextTuple(tuple);
            }

    imageData->GetPointData()->AddArray(grad);
    imageData->Update();


    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName("gradient.vti");
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();
}

void NarrowBand::computeInitialSurface() {
    std::clock_t t0 = clock();
    PriorityQueue boundary_pq;

    for(int k = 1; k < _vDim-1;k++)
        for(int j = 1; j < _vDim-1;j++)
            for(int i = 1; i < _vDim-1;i++) {
                int ijk = fromIJK2int(i,j,k,_vDim);
                _U[ijk] = -1;
                if(i==1 || i==(_vDim-2) || j==1 || j==(_vDim-2) || k==1 || k==(_vDim-2) ) {
                    boundary_pq.push(std::pair<int,float>(ijk,_D[ijk]));
                    _U[ijk] = 0;
                }
            }

    //    PriorityQueue dummy(boundary_pq);
    //    int last = 0;
    //    while(!dummy.empty())
    //    {
    //        last = dummy.top().first;
    //        dummy.pop();
    //    }
    //    writeUtoFile();

    bool stop = false;
    while(stop == false) {
        int grid[3];
        std::pair<int,float> top = boundary_pq.top();
        int current = top.first;
        fromint2IJK(current,_vDim,grid);

        std::vector<std::pair<int,float> > neighbors;
        bool boundary_point = false;
        for(int c = -1; c < 2; c++)
            for(int b = -1; b < 2; b++)
                for(int a = -1; a < 2; a++) {
                    int idx = fromIJK2int(grid[0]+a,grid[1]+b,grid[2]+c,_vDim);
                    if(_U[idx] == -1) {
                        if(_D[idx] <= _D[current])
                            neighbors.push_back(std::pair<int,float>(idx,_D[idx]));
                        else
                            boundary_point = true;
                    }
                }

        if(boundary_point == true) {
            _U[current] = 0;
            boundary_pq.pop();
        }
        else {
            _U[current] = 1;
            boundary_pq.pop();
            for(std::vector<std::pair<int,float> >::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
            {
                std::pair<int,float> temp = *it;
                boundary_pq.push(temp);
                _U[temp.first] = 0;
            }
        }

        //        if(current == last)
        //        {
        //            dummy = boundary_pq;
        //            while(!dummy.empty())
        //            {
        //                last = dummy.top().first;
        //                dummy.pop();
        //            }
        //            _iter++;
        //            writeUtoFile();
        //        }

        if(top.second < 2*_delta)
            stop = true;
    }


    while(!boundary_pq.empty()) {
        _U[boundary_pq.top().first] = 0;
        boundary_pq.pop();
    }

    std::clock_t t1 = clock();
    double elapsed_time1 = double(t1 - t0)/CLOCKS_PER_SEC;
    std::cout << "\tTime to compute Initial Guess: " << elapsed_time1 << "s\n";
}

int NarrowBand::tagCell(int i, int j, int k) {
    int interior=0,exterior=0;
    for(int c = 0; c < 2; c ++)
        for(int b = 0; b < 2; b ++)
            for(int a = 0; a < 2; a ++) {
                int v_ijk = fromIJK2int(i+a,j+b,k+c,_vDim);
                if(_U[v_ijk] < 0)
                    interior++;
                if(_U[v_ijk] > 0)
                    exterior++;
            }
    if(interior > exterior)
        return -1;
    else
        return 1;
}


void NarrowBand::writeUtoFile() {
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(_vDim,_vDim,_vDim);
    imageData->SetOrigin(_min.x(),_min.y(),_min.z());
    imageData->SetSpacing(_delta,_delta,_delta);
    imageData->SetExtent(0,_vDim-1,0,_vDim-1,0,_vDim-1);
#if VTK_MAJOR_VERSION <= 5
    imageData->SetNumberOfScalarComponents(1);
    imageData->SetScalarTypeToFloat();
#else
    imageData->AllocateScalars(VTK_FLOAT, 1);
#endif

    vtkSmartPointer<vtkFloatArray> phi = vtkSmartPointer<vtkFloatArray>::New();
    phi->SetName("zero_level_set");

    for(int idx = 0; idx < _numElem; idx++)
        phi->InsertNextValue(_U[idx]);

    imageData->GetPointData()->AddArray(phi);
    //imageData->GetCellData()->AddArray(phi);
    imageData->Update();

    std::ostringstream file;
    file << "zero_level_set_" << _iter << ".vti";
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(file.str().c_str());
    //std::cerr << "Saving file: " << file.str().c_str() << std::endl;
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();
}


void NarrowBand::evolve(int iter, float deltaT) {
    std::cout << "\t>> Iterations: " << iter << "\n";
    std::cout << "\t>> Time step: " << deltaT << "\n";
    std::cerr << "\t";
    std::clock_t t0 = clock();
    for(_iter = 1; _iter <= iter; _iter++) {
        int if_count = 0, else_count = 0;

        for(int k=0; k < _vDim; k++)
            for(int j=0; j < _vDim; j++)
                for(int i=0; i < _vDim; i++) {
                    int ijk = fromIJK2int(i,j,k,_vDim);
                    if(_vMask[ijk] == 1) {
                        float sol = 0;
                        float DU[3];
                        int ip1jk = fromIJK2int(i+1,j,k,_vDim);
                        int im1jk = fromIJK2int(i-1,j,k,_vDim);
                        int ijp1k = fromIJK2int(i,j+1,k,_vDim);
                        int ijm1k = fromIJK2int(i,j-1,k,_vDim);
                        int ijkp1 = fromIJK2int(i,j,k+1,_vDim);
                        int ijkm1 = fromIJK2int(i,j,k-1,_vDim);

                        if(_vMask[ip1jk] == 1 && _vMask[im1jk] == 1)
                            DU[0] = (_U[ip1jk] - _U[im1jk])/(2*_delta);
                        else
                            DU[0] = std::min((_U[ip1jk]-_U[ijk])/_delta,(_U[ijk]-_U[im1jk])/_delta);
                        if(_vMask[ijp1k] == 1 && _vMask[ijm1k] == 1)
                            DU[1] = (_U[ijp1k] - _U[ijm1k])/(2*_delta);
                        else
                            DU[1] = std::min((_U[ijp1k]-_U[ijk])/_delta,(_U[ijk]-_U[ijm1k])/_delta);
                        if(_vMask[ijkp1] == 1 && _vMask[ijkm1] == 1)
                            DU[2] = (_U[ijkp1] - _U[ijkm1])/(2*_delta);
                        else
                            DU[2] = std::min((_U[ijkp1]-_U[ijk])/_delta,(_U[ijk]-_U[ijkm1])/_delta);

                        float modDU = sqrt(pow(DU[0],2)+pow(DU[1],2)+pow(DU[2],2));
                        float sigma[3][2];
                        if(sqrt(pow(DU[0],2)+pow(DU[2],2)) != 0) {
                            sigma[0][0] = (-DU[2])/sqrt(pow(DU[0],2)+pow(DU[2],2));
                            sigma[1][0] = 0;
                            sigma[2][0] = DU[0]/sqrt(pow(DU[0],2)+pow(DU[2],2));
                            sigma[0][1] = (-DU[0]*DU[1])/(sqrt(pow(DU[0],2)+pow(DU[2],2))*modDU);
                            sigma[1][1] = sqrt(pow(DU[0],2)+pow(DU[2],2))/modDU;
                            sigma[2][1] = (-DU[1]*DU[2])/(sqrt(pow(DU[0],2)+pow(DU[2],2))*modDU);
                        }
                        else {
                            sigma[0][0] = 1;
                            sigma[1][0] = 0;
                            sigma[2][0] = 0;
                            sigma[0][1] = 0;
                            sigma[1][1] = 0;
                            sigma[2][1] = 1;
                        }

                        float d[2][4]={1,-1,1,-1,
                                       1,1,-1,-1};

                        float vertex[3];
                        fromIJK2Vertex(i,j,k,vertex);

                        //                        float xp = vertex[0] + _DD[ijk].x()*deltaT;
                        //                        float yp = vertex[1] + _DD[ijk].y()*deltaT;
                        //                        float zp = vertex[2] + _DD[ijk].z()*deltaT;

                        for(int idx=0; idx < 4; idx++) {

                            float xp = vertex[0] + _DD[ijk].x()*deltaT + sqrt(2*deltaT*_D[ijk])*(sigma[0][0]*d[0][idx]+sigma[0][1]*d[1][idx]);
                            float yp = vertex[1] + _DD[ijk].y()*deltaT + sqrt(2*deltaT*_D[ijk])*(sigma[1][0]*d[0][idx]+sigma[1][1]*d[1][idx]);
                            float zp = vertex[2] + _DD[ijk].z()*deltaT + sqrt(2*deltaT*_D[ijk])*(sigma[2][0]*d[0][idx]+sigma[2][1]*d[1][idx]);

                            int index[3];
                            fromPoint2IJK(xp,yp,zp,index);
                            int id = fromIJK2int(index[0],index[1],index[2],_cDim);
                            if(id >= 0 && id <= _numCell-1)
                            {
                                if_count++;
                                if(_cMask[id] == 1)
                                    sol += interp(xp,yp,zp);
                                else
                                    sol += tagCell(index[0],index[1],index[2]);
                            }
                            else
                            {
                                else_count++;
                                sol += _U[ijk];
                            }
                        }
                        _Unew[ijk] = sol/4;
                    }
                }

        for(int k=0; k < _vDim; k++)
            for(int j=0; j < _vDim; j++)
                for(int i=0; i < _vDim; i++) {
                    int ijk = fromIJK2int(i,j,k,_vDim);
                    if(_vMask[ijk] == 1)
                        _U[ijk] = _Unew[ijk];
                }
        writeUtoFile();

    }
    std::cerr << "\n";
    std::clock_t t1 = clock();
    double elapsed_time1 = double(t1 - t0)/CLOCKS_PER_SEC;
    std::cout << "\tTime to compute Surface Evolution: " << elapsed_time1 << "s\n";
}

float NarrowBand::interp(float xp, float yp, float zp) {

    int grid[3];
    fromPoint2IJK(xp,yp,zp,grid);

    float vertex[3];
    fromIJK2Vertex(grid[0],grid[1],grid[2],vertex);

    int IJK = fromIJK2int(grid[0],grid[1],grid[2],_vDim);
    int Ip1JK = fromIJK2int(grid[0]+1,grid[1],grid[2],_vDim);
    int IJp1K = fromIJK2int(grid[0],grid[1]+1,grid[2],_vDim);
    int Ip1Jp1K = fromIJK2int(grid[0]+1,grid[1]+1,grid[2],_vDim);
    int IJKp1 = fromIJK2int(grid[0],grid[1],grid[2]+1,_vDim);
    int Ip1JKp1 = fromIJK2int(grid[0]+1,grid[1],grid[2]+1,_vDim);
    int IJp1Kp1 = fromIJK2int(grid[0],grid[1]+1,grid[2]+1,_vDim);
    int Ip1Jp1Kp1 = fromIJK2int(grid[0]+1,grid[1]+1,grid[2]+1,_vDim);

    float u00 = (xp - vertex[0])*(1/_delta)*(_U[Ip1JK] - _U[IJK]) + _U[IJK];
    float u01 = (xp - vertex[0])*(1/_delta)*(_U[Ip1Jp1K] - _U[IJp1K]) + _U[IJp1K];
    float u0  = (yp - vertex[1])*(1/_delta)*(u01 - u00) + u00;

    float u10 = (xp - vertex[0])*(1/_delta)*(_U[Ip1JKp1] - _U[IJKp1]) + _U[IJKp1];
    float u11 = (xp - vertex[0])*(1/_delta)*(_U[Ip1Jp1Kp1] - _U[IJp1Kp1]) + _U[IJp1Kp1];
    float u1  = (yp - vertex[1])*(1/_delta)*(u11 - u10) + u10;

    return (zp - vertex[2])*(1/_delta)*(u1 - u0) + u0;

}

double NarrowBand::computeCloudResolution()
{
    double res = 0.0;
    int n_points = 0;

    std::vector<int> indices (2);
    std::vector<float> sqr_distances (2);
    pcl::KdTreeFLANN<PointType> tree;
    tree.setInputCloud (_cloud);

    for (size_t i = 0; i < _cloud->size (); ++i)
    {
        if (! pcl_isfinite ((*_cloud)[i].x))
        {
            continue;
        }
        //Considering the second neighbor since the first is the point itself.
        tree.nearestKSearch (i, 2, indices, sqr_distances);
        res += sqrt (sqr_distances[1]);
        ++n_points;

    }
    if (n_points != 0)
    {
        res /= n_points;
    }
    return res;
}

float NarrowBand::solve(std::vector<float> &a)
{
    float temp = a[0] + _delta;
    if (temp <= a[1])
    {
        return temp;
    }
    else
    {
        temp = (a[0] + a[1] + sqrt(2*pow(_delta,2)-pow(a[0]-a[1],2)))/2;
        if (temp <= a[2])
        {
            return temp;
        }
        else
            return ( a[0]+a[1]+a[2] + sqrt(pow(a[0]+a[1]+a[2],2) - 3*(pow(a[0],2)+pow(a[1],2)+pow(a[2],2)-pow(_delta,2))) )/3;
    }
}

int NarrowBand::fromIJK2int(int i, int j, int k, int dim)
{
    return (i + dim*j + dim*dim*k);
}

void NarrowBand::fromIJK2Center(int i, int j, int k, float out[])
{
    out[0] = _min.x() + i*_delta + _delta/2;
    out[1] = _min.y() + j*_delta + _delta/2;
    out[2] = _min.z() + k*_delta + _delta/2;
}

void NarrowBand::fromIJK2Vertex(int i, int j, int k, float out[])
{
    out[0] = _min.x() + i*_delta;
    out[1] = _min.y() + j*_delta;
    out[2] = _min.z() + k*_delta;
}

void NarrowBand::fromint2IJK(int in, int dim, int out[])
{
    div_t divresult = div(in,dim);
    out[0] = divresult.rem;
    divresult = div(divresult.quot,dim);
    out[1] = divresult.rem;
    divresult = div(divresult.quot,dim);
    out[2] = divresult.rem;
}

void NarrowBand::fromCenter2Cell(float x, float y, float z, int out[])
{
    out[0] = (x - _delta/2 - _min.x())/_delta;
    out[1] = (y - _delta/2 - _min.y())/_delta;
    out[2] = (z - _delta/2 - _min.z())/_delta;
}

void NarrowBand::fromPoint2IJK(float x, float y, float z, int out[])
{
    //    float xc = x - _min.x(), yc = y - _min.y(), zc = z - _min.z();
    //    float xt = xc - fmod(xc,_delta), yt = yc - fmod(yc,_delta), zt = zc - fmod(zc,_delta);

    //    out[0] = xt/_delta;
    //    out[1] = yt/_delta;
    //    out[2] = zt/_delta;
    out[0] = (int)floor((x-_min.x())/_delta);
    out[1] = (int)floor((y-_min.y())/_delta);
    out[2] = (int)floor((z-_min.z())/_delta);
}

void NarrowBand::fromVertex2IJK(float x, float y, float z, int out[])
{
    out[0] = (x - _min.x())/_delta;
    out[1] = (y - _min.y())/_delta;
    out[2] = (z - _min.z())/_delta;
}

void manageDirectories(std::string directory)
{
    struct stat info;
    if( stat( directory.c_str(), &info ) != 0 )
    {
        mkdir(directory.c_str(),0700);
    }
    else if( info.st_mode & S_IFDIR )
    {
        if(isDirectoryEmpty(directory.c_str()) == 0)
        {
            DIR *theFolder = opendir(directory.c_str());
            struct dirent *next_file;
            char filepath[256];

            while ( (next_file = readdir(theFolder)) != NULL )
            {
                if (0==strcmp(next_file->d_name, ".") || 0==strcmp(next_file->d_name, "..")) { continue; }
                char cwd[1024];
                getcwd(cwd, sizeof(cwd));
                sprintf(filepath, "%s/%s", directory.c_str(), next_file->d_name);
                remove(filepath);
            }
            closedir(theFolder);
        }
    }
    else
        printf( "%s is no directory\n", directory.c_str() );

    chdir(directory.c_str());

}

int isDirectoryEmpty(const char *dirname)
{
    int n = 0;
    struct dirent *d;
    DIR *dir = opendir(dirname);
    if (dir == NULL) //Not a directory or doesn't exist
        return 1;
    while ((d = readdir(dir)) != NULL) {
        if(++n > 2)
            break;
    }
    closedir(dir);
    if (n <= 2) //Directory Empty
        return 1;
    else
        return 0;
}
