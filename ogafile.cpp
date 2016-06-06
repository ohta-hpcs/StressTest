#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <time.h>

/*
const size_t	gmin	     = 20;
const size_t	gmax	     = 30;
const size_t	lmin	     = 5;
const size_t	lmax	     = 20;
const size_t	mod_nmax     = 10000;
const size_t	mod_minsize  = 512;
const size_t	mod_maxsize  = 1536;
const size_t	mod_dummynum = 100000;
*/
size_t	gmin;
size_t	gmax;
size_t	lmin;
size_t	lmax;
size_t	mod_nmax;
size_t	mod_minsize;
size_t	mod_maxsize;
size_t	mod_dummynum;

extern "C" size_t write_dummymod(const std::string filename,size_t gmin, size_t gmax, size_t lmin, size_t lmax, size_t mod_nmax, size_t mod_minsize, size_t mod_maxsize,size_t mod_dummynum);

extern "C" int ogafile(size_t gmin, size_t gmax, size_t lmin, size_t lmax, size_t mod_nmax, size_t mod_minsize, size_t mod_maxsize,size_t mod_dummynum){
	// Initialization for rand()
	long int tmptime;
	time(&tmptime);
	srandom(static_cast<unsigned int>(tmptime));

	// create filename
	std::ostringstream os;
	size_t i = 0;
	while(i < mod_dummynum){
		os.str("");			// buffer clear
		/*
		if(argc > 1){
			std::string tmpstr = argv[1];
			os << argv[1] << "/bi_" << i << ".txt";
		} else
		*/
			os << "bi_" << i << ".txt";

		size_t fsize = write_dummymod(os.str(),gmin, gmax, lmin, lmax, mod_nmax, mod_minsize, mod_maxsize,mod_dummynum);
#ifdef DEBUG
		std::cerr << "Write file : " << os.str() << "\t";
		std::cerr << "Write size : " << fsize << std::endl;
#endif	// DEBUG
		++i;
	}

	return 0;
}

extern "C" size_t write_dummymod(const std::string filename, size_t gmin, size_t gmax, size_t lmin, size_t lmax, size_t mod_nmax, size_t mod_minsize, size_t mod_maxsize,size_t mod_dummynum){
	size_t filesize = 0;

	// Create Dummy Module File
	std::ostringstream os2;
	os2.str("");			// Initialization

	while(true){
		os2.str("");
		size_t gene = (size_t)random() % gmax;
		size_t line = (size_t)random() % lmax;
		if ( gene < gmin) gene = gmin;
		if ( line < lmin) line = lmin;
    	
		// Header
		os2 << "Sample\tName\t";
		for(size_t i = 0; i < gene; ++i)
			os2 << "Rn." << rand() % mod_nmax << ".1\t";
			os2 << "\n";
    
		int dummydata = (int)(random() % 5);
		for(size_t j = 0; j < line; ++j){
			os2 << "Dummy\tNA\t";
			for(size_t i = 0; i < gene; ++i){
				if ( (i % 2) == 0 ) os2 << 0 - dummydata;
				else                os2 << dummydata;
				os2 << "\t";
			}	// End of gene loop
			os2 << "\n";
		}	// End of line loop

		// File Size check
		if( os2.str().size() < mod_minsize || os2.str().size() > mod_maxsize ){
#ifdef DEBUG
			//    std::cerr << "fsize over : " << os2.str().size() << "\n";
			//    std::cerr << os2.str() << std::endl;
			//    exit(1);
#endif	// DEBUG
			continue;
		} else
			break;
	}	// End of While loop
	filesize = os2.str().size() + 1;

	// Create Dummy data
	std::ofstream fout;
	fout.open(filename.c_str(), std::ios::out | std::ios::trunc);
	fout << os2.str();
	fout.close();

	return filesize;
}

