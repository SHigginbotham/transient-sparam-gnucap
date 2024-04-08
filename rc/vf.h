/*
    RECURSIVE CONVOLUTION COMPANION MODEL FOR GNUCAP
    USES VECTOR FITTING ALGORITHM TO FIT S-PARAMETER DATA
    AND ALLOW TRANSIENT SIMULATION OF 1-PORT, LTI S-PARAMETER BLOCKS

    THIS IS THE VECTOR FITTING ALGORITHM HEADER FILE
    IT EMPLOYS THE LAPACK LINEAR ALGEBRA LIBRARY

    Created by:     Sean Higginbotham for M.A.I project
                    Supervisor: Dr. Justin King
                    Department of Electronic and Electrical Engineering,
                    Trinity College Dublin, Ireland, September 2023 - April 2024
*/

#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <filesystem>
#include <chrono>
#include <ctime>

// LAPACK
#include <lapacke.h>
#include <cblas.h>

// OUTPUT OF THIS FUNCTION SHOULD BE 0 IF FAILED, AND 1 IF SUCCESSFUL

using namespace std::complex_literals;

std::complex<double> one(1.0, 0);
std::complex<double> imag = 1i;

bool do_vector_fitting(std::vector<std::complex<double>>& p, std::vector<std::complex<double>>& r, double& rem, int nump, int numi, bool vflog)
{
    auto start = std::chrono::system_clock::now();          // start time of VF run

    // set up log file for VF
    if (vflog)
    {
        try { std::filesystem::remove("vf_log.txt"); }
        catch(...) {}
    }
    std::ofstream log_file;
    if (vflog) log_file.open("vf_log.txt", std::ios::app | std::ios::out);
    std::time_t in_time = std::chrono::system_clock::to_time_t(start);
    if (log_file) log_file << "VECTOR FITTING LOG FILE FOR TR RUN AT : " << std::ctime(&in_time) << "\t====================================\n";

    /////////////////
    // Definitions //
    /////////////////

    // the frequencies at which the data points occur (in Hz)
    std::vector<double> freqs_hz;
    // the data points, are a complex number (cause is freq domain)
    std::vector<std::complex<double>> sp_points;
    // name of the S-param data file
    std::string sp_file = "s_param_data.txt";

    //////////////////////
    // Get S-param data //
    //////////////////////

    in_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    if (log_file) log_file << std::ctime(&in_time) << "\tReading in S-param data...\n";
    std::ifstream s_param_file(sp_file);
    if (s_param_file.is_open())
    {
        std::string row;
        while (getline(s_param_file, row))
        {
            if (row[0] == '!') ;
            else
            {
                char c;
                int8_t i = 0;
                size_t p;
                double real = 0., imag = 0.;
                while (c = row[0])
                {
                    if (isspace(c) || c==',') row.erase(0,1);
                    else
                    {
                        i++;    // 1 = freqs, 2 == real data, 3 == imag data
                        p = row.find_first_of("\t\v\f\n\r ");
                        std::string t = row.substr(0, p);
                        if (i == 1) freqs_hz.push_back(atof(t.c_str()));
                        else if (i == 2) real = atof(t.c_str());
                        else if (i == 3) imag = atof(t.c_str());
                        row.erase(0,t.length());
                    }
                }
                std::complex<double> z_temp(real, imag);
                sp_points.push_back(z_temp);
            }
        }
        s_param_file.close();
    }
    else
    {
        if (log_file) log_file << "\tS-param data not available! Make sure there is a file named 's_param_data.txt' in the same directory...\n\tExiting run...\n";
        if (log_file) log_file.close();
        return 0;
    }

    //////////////////////////
    // Vector Fitting setup //
    //////////////////////////

    int iters = numi;     // # of iterations of the fitting alg
    int num_poles = nump;  // # of poles to use in the fit

    // transform to be in terms of omega (2 * pi * f) and
    // normalise to range of 0 to 1
    std::vector<double> freqs_hz_omega;
    for (int j = 0; j < freqs_hz.size(); j++)
    {   
        freqs_hz_omega.push_back(freqs_hz[j] * 2 * M_PI);
    }
    double max_freq = *(std::max_element(freqs_hz_omega.begin(), freqs_hz_omega.end()));
    for (int j = 0; j < freqs_hz_omega.size(); j++)
    {
        freqs_hz_omega[j] = freqs_hz_omega[j] / max_freq;
    }

    // distribute poles across 2PI range in S-plane as per section 3.2 of paper
    // this acts as our initial pole guess
    double increment = (freqs_hz_omega.back() - freqs_hz_omega[1]) / (num_poles - 1);
    std::vector<double> freq_range;     // frequencies at which (initial) poles should occur (this is DISTINCT from our s-param freqs, which are given by freqs_hz_omega)
    for (int j = 0; j < num_poles; j++)
    {
        freq_range.push_back(2 * M_PI * (freqs_hz_omega[1] + increment * (double)j));
    }
    std::vector<std::complex<double>> poles_guess;   // need complex poles, not just real
    for (auto x : freq_range)
    {
        double real = -0.01 * x;
        double imag = x;
        std::complex<double> temp(real, imag);
        poles_guess.push_back(temp);
    }
    std::vector<std::complex<double>> freqs_comp;   // wants freqs to be COMPLEX NUMBERS with freq on imag part
    for (auto x : freqs_hz_omega)
    {
        std::complex<double> temp(0, x);
        freqs_comp.push_back(temp);
    }

    //////////////////////////////
    // Vector Fitting Algorithm //
    //////////////////////////////

    // the purpose of the algorithm is to get a good fit, and
    // then we can simply extract the poles, residues, and remainder of this fit
    // to use in the RC companion model

    // open file for outputting poles of algorithm, and insert initial guess

    // delete if exists
    if (vflog)
    {
        try { std::filesystem::remove("pole_guess.txt"); }
        catch(...) {}
    }

    std::ofstream pole_write;
    if (vflog) pole_write.open("pole_guess.txt", std::ios::app | std::ios::out);
    if (pole_write)
    {
        pole_write << "iteration, poles in 2pif (note these should be multiplied by max_freq to get full scale values; max_freq = " << max_freq << ")\ninit";
        for (auto x : poles_guess)
        {
            pole_write << ", " << x.real() << '+' << x.imag() << 'i';
        }
        pole_write << '\n';
    }

    std::vector<std::complex<double>> real_poles;   // real poles
    std::vector<std::complex<double>> comp_poles;   // complex poles
    std::vector<std::complex<double>> all_zeroes;   // zeroes computed from sigma function residues
    std::vector<std::complex<double>> keep_zeroes;  // zeroes computed from sigma function residues that are used for next pole iteration

    double remainder;                               // remainder of the VF fit
    std::vector<std::complex<double>> residues;     // residues of the VF fit

    in_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    if (log_file) log_file << std::ctime(&in_time) << "\tRunning Vector Fitting algorithm with " << num_poles << " starting poles for " << iters << " iterations...\n";

    // PERFORM ITERATIONS TO OBTAIN OPTIMUM FIT

    if (vflog)
    {
        try { std::filesystem::remove("residual_log.txt"); }
        catch(...) {}
    }
    std::ofstream resid_log;
    if (vflog) resid_log.open("residual_log.txt", std::ios::app | std::ios::out);
    // std::time_t in_time = std::chrono::system_clock::to_time_t(start);
    if (resid_log) resid_log << "Output of VF residuals over iterations;\t[iteration, #rows A_real, #cols A_real, rank, residual]\n";

    for (int iteration = 0; iteration <= iters; iteration++)
    {
        // here we use our poles guess, the read-in S-param data, and the frequencies they occur at
        // we form the Ax = B linear system, and solve for the residues and remainder using LAPACKE

            // Steps are as follows:

            //  (i)     form A matrix as per appendix A of VF paper
            //  (ii)    B matrix is just s-param data (sp_points) separated into real and imag parts
            //  (iii)   separate out real and imag parts of A (also as per app. A)
            //  (iv)    solve Ax=B with LAPACKE
            //  (v)     extract zeroes (i.e., poles of next iteration and of final fit) from this solution with LAPACKE, these act as poles of the next iteration
            //  (vi)    iterate steps (1) to (5) for specified # iterations or until convergence
            //  (vii)   when iterations are done, run least squares a final time to get the residues and remainder
            //          corresponding to final set of poles
            //  (viii)  return to RC model with the correct poles, residues, and remainder describing the fit 


        // (1) PRELIMINARIES

        // poles are updated on each iteration, so make sure the updated
        // ones remain stable (-ve real part)
        for (int k = 0; k < poles_guess.size(); k++)
        {
            if (poles_guess[k].real() > 0)
            {
                std::complex<double> temp(poles_guess[k].real() * -1.0, poles_guess[k].imag());
                poles_guess[k] = temp;
            }
        }
        // divide into real and complex poles
        for (auto x : poles_guess)
        {
            if (x.imag() == 0) real_poles.push_back(x);
            else comp_poles.push_back(x);
        }

        // here we define sizes and create arrays. This needs to be done each iteration because the number of poles can CHANGE
        // with each iteration, so we cannot use a fixed set up
        // NOTE
        // s-param data is given by sp_points
        // (scaled to 2PIf) frequencies for sp_points is given by freqs_comp
        // poles are given by poles_guess

        // FREQUENTLY USED SIZES AND THINGS
        // --------------------------------
        std::size_t num_freqs = freqs_comp.size();                      // # S-param frequencies
        std::size_t num_rp = real_poles.size();                         // # real poles
        std::size_t num_cp = comp_poles.size();                         // # complex poles
        std::size_t num_sp = sp_points.size();                          // # S-param points
        std::size_t num_poles = poles_guess.size();
        int pole_iter = 0;
        double residual = 0;

        // VARIABLE ARRAY SIZES
        // --------------------

        std::size_t a_matrix_r = num_freqs,                 a_matrix_c = 2 * (num_rp + 2 * num_cp) + 1;
        std::size_t b_matrix_r = num_sp * 2 - 1,            b_matrix_c = 1;
        std::size_t a_matrix_real_r = 2 * num_freqs - 1,    a_matrix_real_c = 2 * (num_rp + 2 * num_cp) + 1;

        std::size_t Az_real_r = num_rp,                     Az_real_c = num_rp;
        std::size_t bz_real_r = num_rp,                     bz_real_c = 1;
        std::size_t c_real_r = 1,                           c_real_c = num_rp;
        std::size_t bc_real_r = num_rp,                     bc_real_c = num_rp;
        std::size_t H_real_r = num_rp,                      H_real_c = num_rp;
        std::size_t real_zeroes_real_d = num_rp,            real_zeroes_imag_d = num_rp;

        std::size_t Az_comp_r = 2 * num_cp,                 Az_comp_c = 2 * num_cp;
        std::size_t bz_comp_r = 2 * num_cp,                 bz_comp_c = 1;
        std::size_t c_comp_r = 1,                           c_comp_c = 2 * num_cp;
        std::size_t bc_comp_r = 2 * num_cp,                 bc_comp_c = 2 * num_cp;
        std::size_t H_comp_r = 2 * num_cp,                  H_comp_c = 2 * num_cp;
        std::size_t comp_zeroes_real_d = 2 * num_cp,        comp_zeroes_imag_d = 2 * num_cp;

        // for least squares solution (dgelss)
        lapack_int  m_soln = a_matrix_real_r,                           // # rows of A
                    n_soln = a_matrix_real_c,                           // # columns of A
                    nrhs_soln = 1,                                      // # columns of B
                    lda_soln = n_soln,                                  // leading dimension of A (== # of rows for COL-MAJOR, == # cols for ROW-MAJOR)
                    ldb_soln = b_matrix_c,                              // leading dimension of B (== # of rows for COL-MAJOR, == # rows for ROW-MAJOR)
                    info_soln;                                          // returns '0' if is successful
        double singval[n_soln];
        int jpvt[n_soln];
        int rank;

        // for real zeroes matrix multiplication (dgemm)
        lapack_int  m_mult_r = bz_real_r,                               // # rows of bz_real (column vector)
                    n_mult_r = c_real_c,                                // # columns of c_real (row vector)
                    k_mult_r = bz_real_c,                               // # cols bz_real == # rows c_real
                    lda_mult_r = m_mult_r,                              
                    ldb_mult_r = k_mult_r,
                    ldc_mult_r = num_rp;                                // 'first dimension' of output matrix

        // for real zeroes eigenvalues (dgeev)
        lapack_int  info_eig_r,
                    n_eig_r = num_rp,                                   // order of matrix
                    lda_eig_r = n_eig_r;

        // for comp zeroes matrix multiplication (dgemm)
        lapack_int  m_mult_c = bz_comp_r,                               // # rows of bz_comp (column vector)
                    n_mult_c = c_comp_c,                                // # columns of c_comp (row vector)
                    k_mult_c = bz_comp_c,                               // # cols bz_comp == # rows c_comp
                    lda_mult_c = m_mult_c,
                    ldb_mult_c = k_mult_c,
                    ldc_mult_c = 2 * num_cp;                            // first dimension of output matrix

        // for comp zeroes eigenvalues (dgeev)                          
        lapack_int  info_eig_c,
                    n_eig_c = 2 * num_cp,                               // order of matrix (eigenvalues only defined for square matrices, so order is equivalent to #rows == #cols of H matrix)
                    lda_eig_c = n_eig_c;                                // leading dimension of H matrix (== max(1, N))


        // VARIABLE ARRAY DECLARATIONS
        // note they're defined as just 1D array so that we don't need to use pointers to pointers
        // ---------------------------

        std::complex<double> *a_matrix =        new std::complex<double>[a_matrix_r * a_matrix_c];
        double *b_matrix =                      new double[b_matrix_r * b_matrix_c];
        double *b_matrix_final =                new double[(b_matrix_r + 1) * b_matrix_c];    // for final extraction; we don't ignore DC frequency here so add + 1;
        double *a_matrix_real =                 new double[a_matrix_real_r * a_matrix_real_c];

        double *Az_real =                       new double[Az_real_r * Az_real_c];      // should be zero initialised
        double *bz_real =                       new double[bz_real_r * bz_real_c];
        double *c_real =                        new double[c_real_r * c_real_c];
        double *bc_real =                       new double[bc_real_r * bc_real_c];
        double *H_real =                        new double[H_real_r * H_real_c];
        double *real_zeroes_real =              new double[real_zeroes_real_d];
        double *real_zeroes_imag =              new double[real_zeroes_imag_d];

        double *Az_comp =                       new double[Az_comp_r * Az_comp_c];      // should be zero initialised
        double *bz_comp =                       new double[bz_comp_r * bz_comp_c];      // should be zero initialised
        double *c_comp =                        new double[c_comp_r * c_comp_c];
        double *bc_comp =                       new double[bc_comp_r * bc_comp_c];
        double *H_comp =                        new double[H_comp_r * H_comp_c];
        double *comp_zeroes_real =              new double[comp_zeroes_real_d];
        double *comp_zeroes_imag =              new double[comp_zeroes_imag_d];

        // ZERO-INITIALISE APPROPRIATE ARRAYS
        for (int j = 0; j < Az_real_r * Az_real_c; j++) Az_real[j] = 0;
        for (int j = 0; j < Az_comp_r * Az_comp_c; j++) Az_comp[j] = 0;
        for (int j = 0; j < bz_comp_r * bz_comp_c; j++) bz_comp[j] = 0;

        // (2) PERFORM VECTOR FITTING ITSELF

        in_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        if (log_file) log_file << std::ctime(&in_time) << "\tIteration " << iteration << '\t' << num_rp << '\t' << num_cp << '\n'; // FOR DEBUGGING

        // FILL IN ARRAYS
        // --------------

        // (i) a_matrix
        //  A matrix will be P X 2N + 1
            //  P rows corresponding to frequencies at which s-params were collected (so same size as B matrix)
            //  N columns each side of the '1' term corresponding to specified poles
            //
            //    -left half of A matrix is of form:
            //        1/(s-p) for real poles
            //        1/(s-p) + 1/(s-p'), j/(s-p) - j/(s-p') FOR EACH complex pole (conjugate pair; i.e., there are TWO entries for each complex pole)
            //        All of the first entries come first, then all of the second entries
            //    -right half of A matrix is the same form, but just scaled by the s-param data in the numerator of the partial fractions
            //    -in the middle is a '1' corresponding to the remainder term; note we leave out the proportional term cause supposedly
            //     this is 0 for S-parameters...
        for (int j = 0; j < num_freqs; j++)                                         // rows
        {
            pole_iter = 0;
            for (int k = 0; k < 2 * (num_rp + 2 * num_cp) + 1; k++)       // cols
            {
                if (k >= 0 && k < num_rp)                                            // LHS real entries
                {
                    a_matrix[j*a_matrix_c + k] = one / (freqs_comp[j] - real_poles[pole_iter]);

                    // check for issues
                    if (std::isinf(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].imag()) || std::isinf(a_matrix[j*a_matrix_c + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in LHS real\nExiting run..." << std::endl;
                        if (log_file) log_file.close();
                        if (pole_write) pole_write.close();
                        return 0; 
                    }
                    pole_iter++;
                }
                else if (k >= num_rp && k < num_rp + num_cp)   // LHS comp poles, corresponding to real part of residue for the pair
                {
                    a_matrix[j*a_matrix_c + k] = one / (freqs_comp[j] - comp_poles[pole_iter - num_rp]) + one / (freqs_comp[j] - std::conj(comp_poles[pole_iter - num_rp]));

                    // check for issues
                    if (std::isinf(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].imag()) || std::isinf(a_matrix[j*a_matrix_c + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in LHS comp 1 with pole_iter of " << pole_iter << std::endl;
                        if (log_file) log_file.close();
                        if (pole_write) pole_write.close();
                        return 0;  
                    }
                    pole_iter++;
                }
                else if (k >= num_rp + num_cp && k < num_rp + 2 * num_cp) // LHS comp poles, corresponding to imag part of residue for the pair
                {
                    a_matrix[j*a_matrix_c + k] = imag * ((one / (freqs_comp[j] - comp_poles[pole_iter - num_rp - num_cp])) - (one / (freqs_comp[j] - std::conj(comp_poles[pole_iter - num_rp - num_cp]))));

                    // check for issues
                    if (std::isinf(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].imag()) || std::isinf(a_matrix[j*a_matrix_c + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in LHS comp 2 with pole_iter of" << pole_iter << std::endl; 
                        if (log_file) log_file.close();
                        if (pole_write) pole_write.close();
                        return 0;  
                    }
                    pole_iter++;
                }
                else if (k == num_rp + 2 * num_cp)                        // insert '1' for remainder
                {
                    a_matrix[j*a_matrix_c + k] = one;

                    // check for issues
                    if (std::isinf(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].imag()) || std::isinf(a_matrix[j*a_matrix_c + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in '1' column" << std::endl;
                        if (log_file) log_file.close();
                        if (pole_write) pole_write.close();
                        return 0;  
                    }
                    pole_iter = 0;                                                              // reset pole_iter for RHS
                }
                else if (k > num_rp + 2 * num_cp && k <= num_rp + 2 * num_cp + num_rp )       // RHS real poles
                {
                    a_matrix[j*a_matrix_c + k] = -sp_points[j] / (freqs_comp[j] - real_poles[pole_iter]);

                    // check for issues
                    if (std::isinf(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].imag()) || std::isinf(a_matrix[j*a_matrix_c + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in RHS real" << std::endl; 
                        if (log_file) log_file.close();
                        if (pole_write) pole_write.close();
                        return 0;  
                    }
                    pole_iter++;
                }
                else if (k > num_rp + 2 * num_cp + num_rp && k <= num_rp + 2 * num_cp + num_rp + num_cp)   // RHS comp poles, corresponding to real part of residue for the pair
                {
                    a_matrix[j*a_matrix_c + k] = -sp_points[j] / (freqs_comp[j] - comp_poles[pole_iter - num_rp]) + (-sp_points[j]) / (freqs_comp[j] - std::conj(comp_poles[pole_iter - num_rp]));

                    // check for issues
                    if (std::isinf(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].imag()) || std::isinf(a_matrix[j*a_matrix_c + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in RHS comp 1 with pole_iter of " << pole_iter << std::endl;
                        if (log_file) log_file.close();
                        if (pole_write) pole_write.close();
                        return 0;  
                    }
                    pole_iter++;
                }  
                else if (k > 2 * num_rp + 3 * num_cp && k <= num_rp + 2 * num_cp + num_rp + 2 * num_cp)     // RHS comp poles, corresponding to imag part of residue for the pair
                {
                    a_matrix[j*a_matrix_c + k] = -(imag * sp_points[j]) / (freqs_comp[j] - comp_poles[pole_iter - num_rp - num_cp]) - (imag * (-sp_points[j])) / (freqs_comp[j] - std::conj(comp_poles[pole_iter - num_rp - num_cp]));

                    // check for issues
                    if (std::isinf(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].real()) || std::isnan(a_matrix[j*a_matrix_c + k].imag()) || std::isinf(a_matrix[j*a_matrix_c + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in RHS comp 2 with pole_iter of " << pole_iter << std::endl;
                        if (log_file) log_file.close();
                        if (pole_write) pole_write.close();
                        return 0;  
                    }
                    pole_iter++;
                }
            }
        }

        // (ii) b_matrix
        for (int k = 0; k < num_sp; k++)
        {
            b_matrix[k*b_matrix_c + 0] = sp_points[k].real();
            b_matrix_final[k*b_matrix_c + 0] = sp_points[k].real();
            b_matrix_final[(num_sp + k)*b_matrix_c + 0] = sp_points[k].imag();    // for final don't skip DC
            
            // check for issues
            if (std::isinf(b_matrix[k*b_matrix_c + 0]) || std::isnan(b_matrix[k*b_matrix_c + 0]))
            {
                if (log_file) log_file << "inf or NaN in index " << k << " in B matrix\nExiting run...";
                if (log_file) log_file.close();
                if (pole_write) pole_write.close();
                return 0;
            }
        }
        for (int k = 0; k < (num_sp - 1); k++)
        {
            b_matrix[(num_sp + k)*b_matrix_c + 0] = sp_points[1+k].imag();  // skip DC component for imag part

            // check for issues
            if (std::isinf(b_matrix[(num_sp + k)*b_matrix_c + 0]) || std::isnan(b_matrix[(num_sp + k)*b_matrix_c + 0]))
            {
                if (log_file) log_file << "inf or NaN in index " << k << " in B matrix\nExiting run...";
                if (log_file) log_file.close();
                if (pole_write) pole_write.close();
                return 0;
            }
        }

        // (iii) a_matrix_real
        // system matrix A decomposed into real and imag parts, minus 1 in rows cause we skip DC component for imag half
        for (int j = 0; j < num_freqs; j++)                                         // rows
        {
            for (int k = 0; k < 2 * (num_rp + 2 * num_cp) + 1; k++)       // cols
            {
                a_matrix_real[j*a_matrix_real_c + k] = a_matrix[j*a_matrix_c + k].real();
            }
        }
        for (int j = 0; j < num_freqs - 1; j++)                                         // rows (-1 cause skipping DC component)
        {
            for (int k = 0; k < 2 * (num_rp + 2 * num_cp) + 1; k++)       // cols
            {
                a_matrix_real[(j+num_freqs)*a_matrix_real_c + k] = a_matrix[(1+j)*a_matrix_c + k].imag();            // skip DC component for rows
            }
        }

        // (iv) solve Ax = B using LAPACKE (should probably move this and all above into first iteration loop?)

        // this function (dgelss) solves a linear system where system matrix (A) may be rank deficient
        // solution matrix will have #rows = #cols of A and #cols = #cols of B, and be stored in 'b_matrix'
        // residual of the solution is the sum of squares of the n + 1 : m elements in b_matrix
        // but this is supposedly only 'valid' if m_soln > n_soln && rank == n_soln
        info_soln = LAPACKE_dgelss(LAPACK_ROW_MAJOR, m_soln, n_soln, nrhs_soln, a_matrix_real, lda_soln, b_matrix, ldb_soln, singval, -1, &rank);
        
        if (info_soln != 0)
        {
            if (log_file) log_file << "LAPACKE_dgelss failed on iteration << " << iteration << "with return value: " << info_soln << "\nExiting run...\n";
            if (pole_write) pole_write.close();
            if (log_file) log_file.close();
            return 0;
        }

        // THIS GIVES (SUM OF SQUARES) RESIDUALS OF THE VF SYSTEM, CAN USE FOR COMPARING TO ACTUAL DATA IF DESIRED
        // if (m_soln > n_soln && rank == n_soln)
        // {
            for (int j = a_matrix_real_c; j < a_matrix_real_r; j++)
            {
                for (int k = 0; k < b_matrix_c; k++)
                {
                    residual += std::pow(b_matrix[j * b_matrix_c + k], 2);
                }
            }
            // if (resid_log) resid_log << iteration << ", " << m_soln << ", " << n_soln << ", " << rank << ',' << residual << '\n';
        // }
        /*else */if (resid_log) resid_log << iteration << ", " << m_soln << ", " << n_soln << ", " << rank << ',' << residual << '\n';

        // (v) EXTRACT ZEROES FROM RESIDUES OF THE SIGMA FUNCTION

        // because C is 0 indexed, we use the final iteration to get the final fit (step (vii))
        if (iteration < iters) 
        {

            // We construct the H matrix for real poles/residues and solve for eigenvalues (real zeroes of sigma function)
            // to get poles of next iteration

            for (int k = 0; k < num_rp; k++)
            {
                Az_real[k*Az_real_c + k] = real_poles[k].real();        // diagonal matrix of poles
                // check for issues
                if (std::isnan(Az_real[k*Az_real_c + k]) || std::isinf(Az_real[k*Az_real_c + k]))
                {
                    if (log_file) log_file << "inf or Nan at index (" << k << ',' << k << ") of Az_real!\n";
                    if (pole_write) pole_write.close();
                    if (log_file) log_file.close();
                    return 0;
                }
                
                bz_real[k*bz_real_c + 0] = 1;                       // column vector of 1s

                c_real[0*c_real_c + k] = b_matrix[(num_rp + 2 * num_cp + 1 + k)*b_matrix_c + 0];    // column vector of sigma residues
                // check for issues
                if (std::isnan(c_real[k]) || std::isinf(c_real[k]))
                {
                    if (log_file) log_file << "inf or Nan at index " << k << " of c_real!\n";
                    if (pole_write) pole_write.close();
                    if (log_file) log_file.close();
                    return 0;
                }
            }
            if (num_rp > 0)      // some params of cblas_dgemm must be at least 1, so if no real poles exist then this operation won't be valid and should skip it
            {
                // matrix multiplication to obtain bz * c, and then elementwise subtraction to get the H matrix
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m_mult_r, n_mult_r, k_mult_r, 1.0, bz_real, lda_mult_r, c_real, ldb_mult_r, 0, bc_real, ldc_mult_r);
                for (int k = 0; k < num_rp; k++)
                    for (int j = 0; j < num_rp; j++)
                        H_real[k*H_real_c + j] = Az_real[k*Az_real_c + j] - bc_real[j*bc_real_c + k];       // bc_real is stored as transpose of what we want, so swap row and col indices when accessing to maintain correct operations

                // GET EIGENVALUES for real poles; real parts are stored in real_zeroes_real, and imag parts in real_zeroes_imag 
                info_eig_r = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', n_eig_r, H_real, lda_eig_r, real_zeroes_real, real_zeroes_imag, NULL, 1, NULL, 1);
                if (info_eig_r != 0)
                {
                    if (log_file) log_file << "LAPACKE_dgeev failed with real poles!\nExiting run...\n";
                    if (pole_write) pole_write.close();
                    if (log_file) log_file.close();
                    return 0;
                }
            }

            // construct H matrix for complex poles/residues and solve for eigenvalues (complex zeroes of sigma function)

            for (int k = 0; k < num_cp; k++)
            {
                Az_comp[k*Az_comp_c + k] = comp_poles[k].real();                            // top left sub matrix (real coeffs)
                Az_comp[k*Az_comp_c + (num_cp + k)] = comp_poles[k].imag();                 // top right sub matrix (imag coeffs)
                Az_comp[(num_cp + k)*Az_comp_c + k] = -comp_poles[k].imag();                // bottom left sub matrix (-ve of imag coeffs)
                Az_comp[(num_cp + k)*Az_comp_c + (num_cp + k)] = comp_poles[k].real();      // bottom right sub matrix (real coeffs)
                bz_comp[k*bz_comp_c + 0] = 2;               // first half is a bunch of 2s, rest is 0s
                c_comp[0*c_comp_c + k]  = b_matrix[(2 * num_rp + 2 * num_cp + 1 + k)*b_matrix_c + 0];   // 1st half is real part of residues
                c_comp[0*c_comp_c + (num_cp + k)] = b_matrix[(2 * num_rp + 2 * num_cp + 1 + num_cp + k)*b_matrix_c + 0];    // 2nd half is imag part of residues
            }
            if (num_cp > 0)
            {
                // matrix multiplication to obtain bz * c, and then elementwise subtraction to get the H matrix
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m_mult_c, n_mult_c, k_mult_c, 1.0, bz_comp, lda_mult_c, c_comp, ldb_mult_c, 0, bc_comp, ldc_mult_c);
                for (int k = 0; k < 2 * num_cp; k++)
                {
                    for (int j = 0; j < 2 * num_cp; j++)
                    {
                        H_comp[k*H_comp_c + j] = Az_comp[k*Az_comp_c + j] - bc_comp[j*bc_comp_c + k];       // bc_comp is stored as the transpose of what we want, so swap row and col indices when accessing
                    }
                }

                // GET EIGENVALUES for complex poles; real parts are stored in comp_zeroes_real, and imag parts in comp_zeroes_imag 
                info_eig_c = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', n_eig_c, H_comp, lda_eig_c, comp_zeroes_real, comp_zeroes_imag, NULL, 1, NULL, 1);
                if (info_eig_c != 0)
                {
                    if (log_file) log_file << "LAPACKE_dgeev failed with complex poles!\nExiting run...\n";
                    if (pole_write) pole_write.close();
                    if (log_file) log_file.close();
                    return 0;
                }
            }

            // combine real and complex parts of zeroes;
            for (int k = 0; k < num_rp; k++)
            {
                std::complex<double> temp(real_zeroes_real[k], real_zeroes_imag[k]);
                all_zeroes.push_back(temp);
            }
            for (int k = 0; k < 2 * num_cp; k++)
            {
                std::complex<double> temp(comp_zeroes_real[k], comp_zeroes_imag[k]);
                all_zeroes.push_back(temp);
            }

            // only keep zeroes with positive imag part (one of the complex conj pairs) in 'keep zeroes', and also filter out those with very small magnitudes (meaning they don't contribute to the fit)
            for (auto x : all_zeroes)
            {
                if (x.imag() >= 0)
                {
                    if (std::abs(x) > 1E-10 )
                    {
                        keep_zeroes.push_back(x);
                    }
                }
            }

            // SUB IN ZEROES AS POLES OF NEXT ITERATION

            poles_guess.clear();            // make sure vector is empty (probably doesn't make a difference...)
            poles_guess = keep_zeroes;

        }
        // (vii) after done iterating, take out the appropriate stuff
        //       we neeed to just get the fit residues and remainder that correspond to the last set of poles used to solve the Ax=B system (otherwise our fit wouldn't match up!)
        else if (iteration == iters)
        {
            // set up
            // since we only do this final iteration once it is MORE efficient to define everthing down here rather than with
            // the PRELIMINARIES above
            int a_matrix_final_r = num_freqs;
            int a_matrix_final_c = 2 * num_poles + 1;
            int a_matrix_final_real_r = 2 * num_freqs;
            int a_matrix_final_real_c = 2 * num_poles + 1;
            std::complex<double> *a_matrix_final = new std::complex<double>[a_matrix_final_r * a_matrix_final_c];
            double *a_matrix_final_real = new double[a_matrix_final_real_r * a_matrix_final_real_c];

            int m_final = a_matrix_final_real_r;
            int n_final = a_matrix_final_real_c;
            int nrhs_final = b_matrix_c;
            int lda_final = n_final;
            int ldb_final = nrhs_final;

            double singval_final[n_final];
            int jpvt_final[n_final];
            int rank_final;

            std::vector<double> residues_real;
            std::vector<double> residues_imag;

            // create a_matrix for final extraction
            // note that since we only care about the fit residues, we don't need the RHS of the system matrix
            for (int j = 0; j < a_matrix_final_r; j++)     // rows
            {
                for (int k = 0; k < num_poles; k++)     // cols (only need to loop over num poles)
                {
                    a_matrix_final[j * (2 * num_poles+1) + k] = one / (freqs_comp[j] - poles_guess[k]) + one / (freqs_comp[j] - std::conj(poles_guess[k]));
                    a_matrix_final[j * (2 * num_poles+1) + (num_poles+k)] = imag * ( one / (freqs_comp[j] - poles_guess[k]) - one / (freqs_comp[j] - std::conj(poles_guess[k])) );

                    // check for issues
                    if (std::isinf(a_matrix_final[j * (2 * num_poles+1) + k].real()) || std::isnan(a_matrix_final[j * (2 * num_poles+1) + k].real()) || std::isnan(a_matrix_final[j * (2 * num_poles+1) + k].imag()) || std::isinf(a_matrix_final[j * (2 * num_poles+1) + k].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << k << " in a_matrix_final norm!\n Exiting run...\n";
                        if (pole_write) pole_write.close();
                        if (log_file) log_file.close();
                        return 0;
                    }
                    if (std::isinf(a_matrix_final[j * (2 * num_poles+1) + (num_poles+k)].real()) || std::isnan(a_matrix_final[j * (2 * num_poles+1) + (num_poles+k)].real()) || std::isnan(a_matrix_final[j * (2 * num_poles+1) + (num_poles+k)].imag()) || std::isinf(a_matrix_final[j * (2 * num_poles+1) + (num_poles+k)].imag()))
                    {
                        if (log_file) log_file << "inf or NaN in index " << j << ", " << num_poles + k << " in a_matrix_final conj!\n Exiting run...\n";
                        if (pole_write) pole_write.close();
                        if (log_file) log_file.close();
                        return 0;
                    }
                }
                a_matrix_final[j * (a_matrix_final_c) + (a_matrix_final_c - 1)] = one;      // 1 for remainder term
            }
            // turn into real matrix
            for (int j = 0; j < a_matrix_final_r; j++)
            {
                for (int k = 0; k < a_matrix_final_real_c; k++)
                {
                    a_matrix_final_real[j * a_matrix_final_real_c + k] = a_matrix_final[j * a_matrix_final_c + k].real();
                    a_matrix_final_real[(j + a_matrix_final_real_r/2) * a_matrix_final_real_c + k] = a_matrix_final[j  * a_matrix_final_c + k].imag();
                }
            }

            // perform a last least squares
            info_soln = LAPACKE_dgelss(LAPACK_ROW_MAJOR, m_final, n_final, nrhs_final, a_matrix_final_real, lda_final, b_matrix_final, ldb_final, singval_final, -1, &rank);
            if (info_soln != 0)
            {
                if (log_file) log_file << "\nLAPACKE_dgelss for final (iteration " << iteration << ") failed: " << info_soln << "\nExiting run...\n";
                if (pole_write) pole_write.close();
                if (log_file) log_file.close();
                return 0;
            }

            // get residues and remainder

            // solution vector is of dimensions n_final x nrhs_final
            // first num_poles terms should be real parts of residues, second num_poles terms should be imag parts
            // final term should be the remainder

            remainder = b_matrix_final[n_final - 1];
            for (int k = 0; k < num_poles; k++)
            {
                residues_real.push_back(b_matrix_final[k]);
                residues_imag.push_back(b_matrix_final[num_poles + k]);
                std::complex<double> temp(b_matrix_final[k], b_matrix_final[num_poles + k]);
                residues.push_back(temp);
            }

            // poles are kept from final iteration; that is, the final computed zeroes become the poles of our final fit

            // delete stuff...
            delete [] a_matrix_final;
            delete [] a_matrix_final_real;
            delete [] b_matrix_final;   
        }

        // print poles to file
        if (pole_write)
        {
            pole_write << "Iteration " <<  iteration;
            for (auto x : poles_guess)
            {
                pole_write << ", " << x.real() << '+' << x.imag() << 'i';
            }
            pole_write << '\n';
        }

        // clear stuff for next iteration (if there is one)
        real_poles.clear();
        comp_poles.clear();
        all_zeroes.clear();
        keep_zeroes.clear();

        // FREE DYNAMIC MEMORY FOR NEXT LOOP
        delete [] a_matrix;
        delete [] b_matrix;
        delete [] a_matrix_real;
        delete [] Az_real;
        delete [] bz_real;
        delete [] c_real;
        delete [] bc_real;
        delete [] H_real;
        delete [] real_zeroes_real;
        delete [] real_zeroes_imag;
        delete [] Az_comp;
        delete [] bz_comp;
        delete [] c_comp;
        delete [] bc_comp;
        delete [] H_comp;
        delete [] comp_zeroes_real;
        delete [] comp_zeroes_imag;   
    }

    in_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    if (log_file) log_file << std::ctime(&in_time) << "\tFinished Vector Fitting Algorithm! Writing out results...\n";


    // open file for outputting poles of final results (for testing; can include in a flag later on whether to do this or not...)
    if (vflog)
    {
        try{ std::filesystem::remove("results.txt"); } catch(...) {}
    }
    std::ofstream results_write;
    if (vflog) results_write.open("results.txt", std::ios::app | std::ios::out);

    // get final stuff to take back to RC model
    // need to scale poles and residues since the algorithm was normalised to range of 0 to 1 for frequency.
    // poles and residues will be in terms of 2 * PI * f (angular frequency), not Hz
    if (results_write) results_write << "\n\t=========================\tPOLES OF FIT (RESCALED by " << max_freq << ", in 2PIf); # poles: " << poles_guess.size() << "\t=============\n";
    for (int k = 0; k < poles_guess.size(); k++)
    {
        poles_guess[k] = poles_guess[k] * max_freq;
        if (results_write) results_write << poles_guess[k].real() << " + " << poles_guess[k].imag() << "i\n";
    }
    if (results_write) results_write << "\n\t============================\tRESIDUES OF FIT (RESCALED, in 2PIf); # residues: " << residues.size() << "\t==============\n";
    for (int k = 0; k < residues.size(); k++)
    {
        residues[k] = residues[k] * max_freq;
        if (results_write) results_write << residues[k].real() << " + " << residues[k].imag() << "i\n";
    }
    if (results_write) results_write << "\n\t==============================\tREMAINDER OF FIT:\t======================\n" << remainder << '\n';

    // finalise and close files
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> time_to_run = end - start;
    in_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    if (log_file) log_file << std::ctime(&in_time) << "\tdo_vector_fitting() finished, it took "<< time_to_run.count() << " seconds to run!\n\nClosing log file...\n";
    if (results_write) results_write.close();
    if (pole_write) pole_write.close();
    if (log_file) log_file.close();
    if (resid_log) resid_log.close();

    // (viii) return fit to RC model via pass by reference and exit VF
    p = poles_guess;
    r = residues;
    rem = remainder;

    return 1;
}
