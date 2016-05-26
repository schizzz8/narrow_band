#include "NarrowBand.h"


std::string filename;
int prec = 5;
int iter = 0;
float deltaT = 0.001;
int expand = 1;


void parseCommandLine(int argc, char* argv[]);


int main(int argc, char* argv[])
{
    parseCommandLine(argc,argv);
    NarrowBand nb(filename);
    nb.build3Dgrid(prec);
    nb.expandGrid(expand);
    nb.writeGridToFile();
    nb.computeDistanceFunction();
    nb.writeDtoFile();
    nb.computeGradient();
    nb.writeDDtoFile();
    nb.computeInitialSurface();
    nb.writeUtoFile();
    nb.evolve(iter,deltaT);

    return 0;
}

void parseCommandLine(int argc, char *argv[])
{
    std::cout << "\n";
    std::cout << BOLD(FRED("Reconstruct implicit surface from a 3D data set with Level-Set method.\n"));
    std::cout << "\n";

    if(argc == 1)
    {
        std::cout << FBLU("Syntax is: narrow_band input.pcd <options>\n");
        std::cout << "  where options are:\n";
        std::cout << "\t-p n" << "\t= grid precision (default: " << FCYN("10") << ")\n";
        std::cout << "\t-e n" << "\t= grid expansion factor (default: " << FCYN("3") << ")\n";
        std::cout << "\t-iter n" << "\t= number of iterations (default: " << FCYN("1") << ")\n";
        std::cout << "\t-dt n" << "\t= time step (default: " << FCYN("0.01") << ")\n";
        exit(0);
    }

    if(argc > 1)
    {
        filename = argv[1];

        for (int i = 2; i < argc; i++) {
            if(strcmp(argv[i],"-p") == 0)
                prec = atoi(argv[i+1]);
            if(strcmp(argv[i],"-e") == 0)
                expand = atoi(argv[i+1]);
            if(strcmp(argv[i],"-iter") == 0)
                iter = atoi(argv[i+1]);
            if(strcmp(argv[i],"-dt") == 0)
                deltaT = strtod(argv[i+1],NULL);
        }
    }
}
