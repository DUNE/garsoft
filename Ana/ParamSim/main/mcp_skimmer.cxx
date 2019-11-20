#include "MCP_Skimmer.h"
#include <iostream>

#include <TFile.h>
#include <TTree.h>

void ShowHelp()
{
    std::cout << "./mcp_skimmer --infile <inputfile> --outfile <outputfile> --debug <0/1>" << std::endl;
}

int main(int argc, char **argv)
{
    if( argc == 1 || ((argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1]))) || argc < 7 ){
        ShowHelp();
        return 2;
    }

    if( argv[1] != std::string("--infile") || argv[3] != std::string("--outfile") || argv[5] != std::string("--debug") ) {
        ShowHelp();
        return -2;
    }

    // get command line options
    std::string outfile = "";
    std::string infile = "";
    std::string debug = "";
    int p = 0;
    while( p < argc )
    {
        if( argv[p] == std::string("--infile") ){
            infile = argv[p+1];
            p++;
        }
        else if( argv[p] == std::string("--outfile") ){
            outfile = argv[p+1];
            p++;
        }
        else if( argv[p] == std::string("--debug") ){
            debug = argv[p+1];
            p++;
        }
        else{
            p++;
        }
    }

    printf( "Skimming from tree dump: %s\n", infile.c_str() );
    printf( "Output Skimmed file: %s\n", outfile.c_str() );
    printf( "Debug mode: %s\n", debug.c_str() );

    MCP_Skimmer *skimmer = new MCP_Skimmer(infile, outfile);
    if(debug == "1")
    skimmer->SetDebug(true);
    if(not skimmer->BookTFile()) return -1;

    skimmer->SkimMCParticle();

    skimmer->WriteTTree();
    skimmer->CloseTFile();
    printf( "-30-\n" );

    return 0;
}
