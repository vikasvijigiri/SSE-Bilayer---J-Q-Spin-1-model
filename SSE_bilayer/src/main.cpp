#include <chrono>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "../include/SSElattice.hpp"
#include "../include/SSEvariables.hpp"
#include "../include/SSEobservables.hpp"
#include "../include/SSEupdates.hpp"			//Headers
#include "../include/writeresults.hpp"


static ran rann;
#define PRINT_VAR(name) std::cout << #name << ": " << vars.name << std::endl




int main() {
    //****************** Creating all necessary instances ********************

    SSEvariables vars;    
    SSEobservables obs;
    SSEupdates update;
    writeresults write;

    vars.declare_variables();

    PRINT_VAR(lx), PRINT_VAR(ly), PRINT_VAR(JH), PRINT_VAR(JQ), PRINT_VAR(JB), PRINT_VAR(Beta),
        PRINT_VAR(isteps), PRINT_VAR(iter), PRINT_VAR(lambda), PRINT_VAR(time_series);

    vars.lattice_sites();

    update.initialize(vars);
    update.weights();



    //****************** define spins info on lattice ********************

    SSElattice* lattice = new SSElattice[vars.Nm];
    for (int xc = 0; xc < vars.Nm; ++xc)
        lattice[xc].set_S((rann() <= 0.5) ? -1 : 1);  // Start with a random configuration



    auto t1 = std::chrono::high_resolution_clock::now();


    // ********************** Equilibration steps *************************
    for (int i = 0; i < vars.isteps; ++i) {
        update.mcstep(vars, lattice);  // Equilibration to check how much time it takes to complete.
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    double time = (elapsed.count() * 1e-9) / 60.;
    std::cout << "Total Equilibration time ~ " << time << " minutes" << std::endl;

    obs.Initiate_observables(vars);  // Initiate Observables.

    // ********************** Measurement steps ***************************
//    std::ofstream outputFile("SSE_data_time_series.txt", std::ios::app);
//    if (!outputFile.is_open()) {
//        std::cerr << "Error opening the file." << std::endl;
//        return 1;
//    }

//    if (outputFile.tellp() == 0)
//        for (const auto& colHeader : {"enrg", "SMag", "stiffx", "stiffy", "O_N", "R_N", "O_V", "R_V", "O_B", "R_B"})
//            outputFile << std::setw(15) << colHeader;

    for (int i = 0; i < vars.iter; ++i) {
        update.mcstep_measurement(vars, lattice, obs);
//        if (vars.time_series)
//            write.save_time_series_data(outputFile, obs.enrg, obs.mag, obs.stiffx, obs.stiffy, obs.O_N, obs.R_N, 
//                        obs.O_V, obs.R_V, obs.O_B, obs.R_B);
    }


    // *************** output printing/processing *************************

//    outputFile.close();



    obs.observables_processing(vars); 

    const char** headers =  obs.createHeadersArray();
    const double* values =  obs.createValuesArray();
    
    write.save_bin_data("SSE_data_avg.txt", headers, values, 33);



    // *********************** Climax *************************************

    t1 = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t2);
    time = (elapsed.count() * 1e-9) / 60.;
    std::cout << "Total Measurements time ~ " << time << " minutes." << std::endl;

    delete[] lattice;
    return 0;
}

