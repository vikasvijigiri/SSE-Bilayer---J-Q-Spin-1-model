//#include "../include/SSEvariables.hpp"
//#include "../include/writeresults.hpp"

//#include <sstream>
//#include <fstream>
//#include <iostream>
//#include <iomanip>
//#include <cstring>

#include <fstream>
#include <iomanip>
#include <cstring>
#include <iostream>
#include <vector>

#include "../include/writeresults.hpp"


using namespace std;




const char* createHeadersArray();


//template <typename... Args>
//void writeresults::save_time_series_data(std::ofstream& outputFile, const Args&... args, const double[3][4], 
//                                        const double[4][4][4], const double[4][2][2]) {
//    ((outputFile << std::setw(15) << std::scientific << args),... );
//    outputFile << '\n';
//}


//***************************************************************************************
/*
template <typename arr, typename... Args>
void writeresults::save_bin_data(arr first, const Args&... args) {


    // Open the file in append mode
    std::ofstream outputFile("SSE_data_avg.txt", std::ios::app);

    // Check if the file is open
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        //return 0; // return an error code
    }

    // Check if the file is empty (to determine whether to write the header)
    outputFile.seekp(0, std::ios::end);
    bool fileIsEmpty = outputFile.tellp() == 0;


    // Write the header if the file is empty
    if (fileIsEmpty) {
      const char* headers[] = {"enrg", "enrg2", "enrg4", 
                               "SMagL1", "SMagL2", 
                               "SMag_squareL1", "SMag_squareL2", "SMag_squareI1", 
                               "SMag_fourL1",  "SMag_fourL2", "SMag_fourI1", "SMag_fourI2", "SMag_fourI3",
                               "stiffx", "stiffy"};

      // Write the header to the file
      // print Headers
      for (const auto& colHeader : headers) {
          outputFile << std::setw(15) << colHeader;
      }
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
          outputFile << std::setw(15) << std::scientific << "N_real" << i << "_" << j;
        }
      }
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          for (int k = 0; k < 4; ++k) {
            outputFile << std::setw(15) << std::scientific << "V_real" << i << "_" << j << "_" << k;
          }        
        }
      }
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
          for (int k = 0; k < 2; ++k) {
            outputFile << std::setw(15) << std::scientific << "B_real" << i << "_" << j << "_";
          }
        }
      }
      outputFile << '\n';
    }

    // print values
    ((outputFile << std::setw(15) << std::scientific << args),... );

    outputFile << '\n';

    // Close the file
    outputFile.close();
}



void 
writeresults::dumpArrayToFile(SSEvariables sv, const double* arr) {
    // Open the file
    std::ofstream outFile("SSE_szsz_avg.txt", std::ios::app);

    // Check if the file is successfully opened
    if (!outFile) {
        std::cerr << "Error opening file: " << "SSE_szsz_avg.txt" << std::endl;
        return;
    }

    // Check if the file is empty (to determine whether to write the header)
    outFile.seekp(0, std::ios::end);
    bool fileIsEmpty = outFile.tellp() == 0;

    // Write the header if the file is empty
    if (fileIsEmpty) {
      for (int i = 0; i < sv.Ns; ++i) {
          outFile << std::setw(15) << std::scientific << "Sz_0Sz_" << i;
      }
      outFile << '\n';
    }


    for (int i = 0; i < sv.Ns; ++i) {
        // Write the element to the current column
        outFile << std::setw(15) << std::scientific << arr[i];
    }
    outFile << '\n';

    // Close the file
    outFile.close();
}
*/



// writeresults.cpp


//template <typename... Args>
void writeresults::save_bin_data(const char* filename , const char** headers , const double* values, int size_values) {
    std::ofstream outputFile(filename, std::ios::app);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return;  // Return an error code or throw an exception if needed
    }

    // Check if the file is empty (to determine whether to write the header)
    outputFile.seekp(0, std::ios::end);
    bool fileIsEmpty = outputFile.tellp() == 0;

    // Write the header if the file is empty
    if (fileIsEmpty) {
        //const char* headers = createHeadersArray();

        // Write the header to the file
    for (int i = 0; i < size_values; i++) {
            outputFile << std::setw(15) << headers[i];
        }

        outputFile << '\n' << '\n';
    }
    for (int i = 0; i < size_values; i++) {
        outputFile << std::setw(15) << std::scientific << values[i];
    }

    // Print values
    //((outputFile << std::setw(15) << std::scientific << args), ...);

    outputFile << '\n';

    // Close the file
    outputFile.close();
}

//template <typename T>
//void writeresults::writeHeader(std::ostream& file, const T& arg) {
//    file << std::setw(15) << std::scientific << arg << '\n';
//}

//template <typename T, typename... Args>
//void writeresults::writeHeader(std::ostream& file, const T& arg, const Args&... args) {
//    file << std::setw(15) << std::scientific << arg;
//    writeHeader(file, args...);
//}

//*******************************************************************************************


//template void writeresults::save_time_series_data(std::ofstream& , const double& , const double& , const double& , const double& 
//, const double& , const double& , const double& , const double& , const double& , const double& , const double& ,
//, const double& , const double& , const double& , const double& , const double[3][4] , const double[4][4][4] , const double[4][2][2]);

//template void writeresults::save_bin_data(const double& , const double& , const double& , const double& 
//, const double& , const double& , const double& , const double& , const double& , const double& , const double& 
//, const double& , const double& , const double& , const double& , const double[3][4] , const double[4][4][4] , const double[4][2][2]);
