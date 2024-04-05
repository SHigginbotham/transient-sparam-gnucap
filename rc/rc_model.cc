/*
    RECURSIVE CONVOLUTION COMPANION MODEL FOR GNUCAP
    USES VECTOR FITTING ALGORITHM TO FIT S-PARAMETER DATA
    AND ALLOW TRANSIENT SIMULATION OF 1-PORT, LTI S-PARAMETER BLOCKS

    Created by:     Sean Higginbotham for M.A.I project
                    Supervisor: Dr. Justin King
                    Department of Electronic and Electrical Engineering,
                    Trinity College Dublin, Ireland, September 2023 - April 2024
*/

// for RC model
#include "globals.h"
#include "e_storag.h"

// for VF
#include "vf.h"

namespace {

    // ================================================
    // = RECURSIVE CONVOLUTION COMPANION MODEL PLUGIN =
    // ================================================

        // THE FPOLY(s) AND CPOLY FOR RC MODEL AND HOW THEY'RE CONSTRUCTED AS TAYLOR SERIES

        // y = (1/Gc) * i_c + (1/Gc) * i[n]                                 Impements the port voltage constitutive relation v[n], is set in do_tr()
        //      y.f0 = (1/Gc)*i_c
        //      y.f1 = 1/Gc
        //      y.x = i[n] = Gc * involts() - i_c

        // m = -i_c + Gc * v[n]                                             Implements the port current constitutive relation i[n], is set in do_tr()
        //      m.c0 = -i_c
        //      m.c1 = Gc
        //      m.x = v[n] = involts()

        // i = i_c = (2 * h_n) / (sqrt(Zref) * (1 + factor))                Implements the i_c RC model internal current, is set in tr_advance()
        //      i.f0 = 0
        //      i.f1 = 2 / (sqrt(Zref) * (1 + factor))
        //      i.x = h_n

    class RC_MODEL : public STORAGE
    {
        private:
            PARAMETER<int> nump;        // # initial poles to use in VF
            PARAMETER<int> numi;        // # iterations to use in VF
            PARAMETER<bool> vflog;      // flag for if want to print logs during the VF (logs are put into 'vf_log.txt', and VF results are put into 'results.txt', and each pole iteration is put into 'pole_guess.txt'')
            PARAMETER<bool> rclog;      // flag for if want to print logs during the transient sims (logs are put into 'tr_write.txt')
        protected:
            explicit RC_MODEL(const RC_MODEL& p) : STORAGE(p) {}         // copy constructor
        public:
            explicit RC_MODEL() : STORAGE() {}                           // constructor

            // parameter function overrides
            int param_count() const             {return (4 + STORAGE::param_count());}  // # params for this device
            bool param_is_printable(int) const;
            std::string param_name(int) const;
            std::string param_name(int, int) const;
            std::string param_value(int) const;
            void set_param_by_index(int i, std::string&, int offset);
        protected:
            char id_letter() const              {return 'p';}
            std::string value_name() const      {return "Zref";}           // pass in Zref
            std::string dev_type() const        {return "RCmodel";}
            int max_nodes() const               {return 2;}
            int min_nodes() const               {return 2;}
            int matrix_nodes() const            {return 2;}
            int net_nodes() const               {return 2;}
            bool has_iv_probe() const           {return 2;}

            // device will be passive cause in this case we're only considering
            // linear, passive (R,L,C) elements in the S-param block

            void tr_iwant_matrix()              {tr_iwant_matrix_passive();}
            void ac_iwant_matrix()              {ac_iwant_matrix_passive();}
            void tr_load()                      {tr_load_passive();}
            void tr_unload()                    {tr_unload_passive();}

            double tr_involts() const           {return tr_outvolts();}
            double tr_involts_limited() const   {return tr_outvolts_limited();}
            COMPLEX ac_involts() const           {return ac_outvolts();}

            // functions which explicitly implement the RC companion model
            bool do_tr();
            void tr_advance();      // update i_c[n] here using _time, _dt, etc...
            void tr_begin();        // for initialisations
            void precalc_first();   // For doing vector fitting and evaluating parameters
            void tr_accept();       // do storage of history terms here 

            // not currently implemented, see how it's done in d_cap.cc
            // can seemingly add ability to probe specific internal states of the device
            // during sims; could use to return i_c[n], etc...
            // double tr_probe_num(const std::string& x) const;

            std::string port_name(int i) const
            {
                assert(i >= 0); assert(i < 2);
                static std::string names[] = {"p", "n"};
                return names[i];
            }

            CARD* clone() const                 {return new RC_MODEL(*this);}

            // FOR STORING DATA FROM VF ALGORITHM

            std::vector<std::complex<double>> vf_poles;
            std::vector<std::complex<double>> vf_residues;
            double vf_remainder;

            // default parameter values (not including Zref)
            const int nump_default = 35;
            const int numi_default = 450;
            const bool vflog_default = 0;
            const bool rclog_default = 0;

            // actual values evaluated
            int num_p;
            int num_i;
            bool vf_log;
            bool rc_log;

            // RC MODEL RELATED PARAMS

            double Zref = 50;                               // Impedance the block is referenced to (Note this is NOT the characteristic impedance of the port TLs)
            double Gc;                                      // parallel conductance of the RC model
            double i_c;                                     // independent current source of the RC model
            double factor;                                  // the (K + 2 * real(sum(rk * lambdak))) term that is repeated in the RC constitutive relations (in Gc and i_c)
            double h_n;                                     // history term ( is actually 2*real(h[n]), where h[n] == hist_sum )
            std::vector<std::complex<double>> dk_store;     // values of the dk[n-1] history terms (must store as vector since will have an entry for each pole)
            double a_store[2] = {0, 0};                     // for storing a[n-1] and a[n-2]
            double a_n;                                     // value of incident wave a[n] at current time step;
            std::vector<std::complex<double>> dk_n;         // value of history term dk[n] at current time step

            std::vector<std::complex<double>> alpha_k, lambda_k, mu_k, nu_k;    // coeficients for 2nd order approx (for time steps > 2)
            std::vector<std::complex<double>> lambda_k_1, mu_k_1;               // 1st order approx (for time steps 1 and 2)
            std::complex<double> res_lambda_sum;                                // dot product of residues and 2nd order lambda
            std::complex<double> res_lambda1_sum;                               // dot product of residues and 1st order lambda
            std::complex<double> hist_sum;                                      // dot product of residues with dk_store and a_store (i.e., the h[n] history term)

            std::ofstream tr_data;  // debugging file for RC model
            bool accepted;          // to prevent tr_accept pushing history terms TWICE in a given time step (which messes up the model, obviously!)

            int time_step;          // current time step (indexed from 0, so time_step = 0 means the 1st time step of the transient run)
    };

    // ======================================================
    // = SET UP PARAMETERS WANNA READ IN WHEN INSTANTIATING =
    // =======================================================

    bool RC_MODEL::param_is_printable(int I) const
    {
        switch (RC_MODEL::param_count() - 1 - I)
        {
            case 0: return nump.has_hard_value();
            case 1: return numi.has_hard_value();
            case 2: return vflog.has_hard_value();
            case 3: return rclog.has_hard_value();
            default: return STORAGE::param_is_printable(I);
        }
    }
    std::string RC_MODEL::param_name(int I) const
    {
        switch (RC_MODEL::param_count() - 1 - I)
        {
            case 0: return "nump";
            case 1: return "numi";
            case 2: return "vflog";
            case 3: return "rclog";
            default: return STORAGE::param_name(I);
        }
    }
    std::string RC_MODEL::param_name(int I, int j) const
    {
        if (j == 0) return param_name(I);
        else if (I >= STORAGE::param_count())
        {
            switch (RC_MODEL::param_count() - 1 - I)
            {
                case 0: return (j==1) ? "nump" : "";
                case 1: return (j==1) ? "numi" : "";
                case 2: return (j==1) ? "vflog" : "";
                case 3: return (j==1) ? "rclog" : "";
                default: return "";
            }
        }
        else return STORAGE::param_name(I, j);
    }
    std::string RC_MODEL::param_value(int I) const
    {
        switch (RC_MODEL::param_count() - 1 - I)
        {
            case 0: return nump.string();
            case 1: return numi.string();
            case 2: return vflog.string();
            case 3: return rclog.string();
            default: return STORAGE::param_value(I);
        }
    }
    void RC_MODEL::set_param_by_index(int I, std::string& Value, int Offset)
    {
        switch (RC_MODEL::param_count() - 1 - I)
        {
            case 0: nump = Value; break;
            case 1: numi = Value; break;
            case 2: vflog = Value; break;
            case 3: rclog = Value; break;
            default: STORAGE::set_param_by_index(I, Value, Offset); break;
        }
    }

    // ===================================
    // = OVVERIDE FUNCTIONS FOR RC MODEL =
    // ===================================

    // is run when performing any command on the device
    void RC_MODEL::precalc_first()
    {
        STORAGE::precalc_first();                       // won't work if don't call base class version first!
        
        // we must first evaluate values of parameters provided with the instance, as otherwise they'd stay as defaults
        num_p = nump.e_val(nump_default, scope());
        num_i = numi.e_val(numi_default, scope());
        vf_log = vflog.e_val(vflog_default, scope());
        rc_log = rclog.e_val(rclog_default, scope());

        // having vector fitting in here seems to prevent issues (i.e., incorrect computatations) during transient sims
        // probably means that the model (device in circuit) needs to be replaced if we change the input s_param_data.txt
        assert(do_vector_fitting(vf_poles, vf_residues, vf_remainder, num_p, num_i, vf_log));
    }

    // initialise relevant values and set up; is run once at beginning of transient sim
    void RC_MODEL::tr_begin()
    {

        STORAGE::tr_begin();

        // open log file for debugging
        try { std::filesystem::remove("tr_write.txt"); } catch(...) {}
        if (rc_log && !tr_data) tr_data.open("tr_write.txt", std::ios::app | std::ios::out);
        if (tr_data) tr_data << "Opening file in tr_begin()...\n";

        // ======================
        // = INITIALIASE VALUES =
        // ======================

        Zref = value();
        i_c = 0;
        time_step = 0;
        accepted = 0;   // for some reason tr_accept() is entered twice, but wanna prevent this!

        _i[0].f0 = 0.;
        _i[0].f1 = 0.;
        _i[0].x = 0.;
        _y[0].f0 = 0.;
        _m0.c0 = 0.;

        // initialise the history terms
        h_n = 0.;
        std::vector<std::complex<double>> zeroes(vf_poles.size(), (0.,0.));
        dk_store = zeroes;  // initialise to zeroes
        a_store[0] = 0.;    // a[n-1]
        a_store[1] = 0.;    // a[n-2]

        if (tr_data) tr_data << "Exiting tr_begin()...\n";
    }

    // is done at start of each time step; here is where the time-step dependant parameters are computed
    // (except for the first time step which is done in do_tr())
    void RC_MODEL::tr_advance()
    {
        // we close tr_data at the end of each time step, so re-open again if needed.
        if (rc_log && !tr_data) tr_data.open("tr_write.txt", std::ios::app | std::ios::out);
        if (tr_data) tr_data << "Entered tr_advance()...\n";
        time_step++;
        STORAGE::tr_advance();

        // clear vectors and things for this iteration
        alpha_k.clear(); lambda_k.clear(); mu_k.clear(); nu_k.clear();
        if (time_step == 1) {lambda_k_1.clear(); mu_k_1.clear();} // 2nd time step (time_step == 1)
        res_lambda_sum = 0. * one;  // set to 0
        if (time_step == 1) res_lambda1_sum = 0. * one;
        hist_sum = 0. * one;    // set to 0

        // make sure VF solution is actually usable!
        assert(vf_poles.size() != 0);
        assert(vf_poles.size() == vf_residues.size());

        // compute alpha, lambda, mu, and nu for this time step
        for (int i = 0; i < vf_poles.size(); i++)           // # poles == # residues
        {
            std::complex<double> x = vf_poles[i];           // pole pk
            std::complex<double> r = vf_residues[i];        // residue rk
            
            std::complex<double>lambda_1, mu_1;
            std::complex<double> alpha = std::exp(x * _dt);
            if (time_step == 1)
            {
                lambda_1 = - (one / x) * ( one + (one - alpha) / (x * _dt) );
                mu_1 = - (one / x) * ( (alpha - one) / (x * _dt) - alpha );
            }
            std::complex<double> lambda = - (one / x) * ((one - alpha)/std::pow((-x * _dt),2) - (3. * one - alpha)/(-2. * x * _dt) + one);
            std::complex<double> mu = - (one / x) * ( (-2. * one) * (one - alpha)/std::pow(-x * _dt,2) + 2. * one /(-x * _dt) - alpha);
            std::complex<double> nu = - (one / x) * ( (one - alpha)/std::pow((- x * _dt),2) - (one + alpha)/(-2. * x * _dt) );

            res_lambda_sum += r * lambda;
            res_lambda1_sum += r * lambda_1;


            // write out log data
            if (time_step == 1)
            {
                tr_data << "\np[" << i+1 << "] = " << x;
                tr_data << "\nr[" << i+1 << "] =" << r;
                tr_data << "\nalpha[" << i+1 << "] = " << alpha;
                tr_data << "\nlambda1[" << i+1 << "] = " << lambda_1;
                tr_data << "\nmu1[" << i+1 << "] = " << mu_1;
                if (tr_data) tr_data << "\nhist_sum[" << time_step + 1 << "] = " << hist_sum << " + ( " << r << " * ( " << alpha << " * " << dk_store[i] << " + " << mu_1 << " * " << a_store[0] << " )) = ";   
                hist_sum += (r * (alpha * dk_store[i] + mu_1 * a_store[0]) );
                if (tr_data) tr_data << hist_sum;
            }
            else
            {   
                tr_data << "\np[" << i+1 << "] = " << x;
                tr_data << "\nr[" << i+1 << "] =" << r;
                tr_data << "\nalpha[" << i+1 << "] = " << alpha;
                tr_data << "\nlambda2[" << i+1 << "] = " << lambda;
                tr_data << "\nmu2[" << i+1 << "] = " << mu;
                tr_data << "\nnu2[" << i+1 << "] = " << nu; 
                if (tr_data) tr_data << "\nhist_sum[" << time_step + 1 << "] = " << hist_sum << " + ( " << r << " * ( " << alpha << " * " << dk_store[i] << " + " << mu << " * " << a_store[0] << " + " << nu << " * " << a_store[1] << " )) = ";   
                hist_sum += (r * (alpha * dk_store[i] + mu * a_store[0] + nu * a_store[1]) );
                if (tr_data) tr_data << hist_sum << "\n";
            }
            alpha_k.push_back(alpha);
            lambda_k.push_back(lambda);
            if (time_step == 1) {lambda_k_1.push_back(lambda_1); mu_k_1.push_back(mu_1);}
            mu_k.push_back(mu);
            nu_k.push_back(nu);
        }

        // ensure we don't end up dividing by 0
        assert(res_lambda_sum.real() != 0);
        if (time_step == 1) assert(res_lambda1_sum.real() != 0);
        assert(vf_remainder != 0);

        factor = (vf_remainder + 2 * (res_lambda_sum.real()));
        double factor_1 = (vf_remainder + 2 * (res_lambda1_sum.real()));    // 1st order for 2nd time step

        // calculate Gc term
        Gc = (1. - factor) / ( Zref * (1. + factor) );                      

        // - calculate h[n] term
        h_n = (2 * (hist_sum.real()));

        // calculate i_c term
        if (time_step == 1) _i[0].f1 = 2 / (std::sqrt(Zref) * (1. + factor_1));
        else _i[0].f1 = 2 / (std::sqrt(Zref) * (1. + factor));
        _i[0].x = h_n;
        i_c = _i[0].f1 * _i[0].x;
        if (tr_data && time_step == 1)
        {
            tr_data << "Time step " << time_step + 1 << '\n';
            tr_data << "\tfactor = " << factor_1 << "\n\tGc = " << Gc << "\n\th_n = " << h_n << "\n\ti_c = " << i_c << '\n'; 
        }
        else if (tr_data)
        {
            tr_data << "Time step " << time_step + 1 << '\n';
            tr_data << "\tfactor = " << factor << "\n\tGc = " << Gc << "\n\th_n = " << h_n << "\n\ti_c = " << i_c << '\n'; 
        }

        if (tr_data) tr_data << "Exiting tr_advance()...\n";
    }

    // do_tr is repeated multiple times till reach satisfied convergence in MNA solver
    bool RC_MODEL::do_tr()
    {

        if (tr_data) tr_data << "Entering do_tr()....\n";
        accepted = 0;       // reset for this time step

        assert(_y[0] == _y[0]);

        // for the first time step, tr_advance() is not entered, so need to do the updating here
        if (time_step == 0)   // first time step (time_step == 0)
        {
            double delta = 0.1e-3;    // assign some arbitrary initial time step of 100ns since initial time step can't be accessed on first step?

            // clear vectors and things for this iteration (note all these will be filled with 1st order stuff as appropriate so no need to discern like in tr_advance())
            alpha_k.clear(); lambda_k.clear(); mu_k.clear(); nu_k.clear();
            res_lambda_sum = 0. * one;  // set to 0
            hist_sum = 0. * one;        // set to 0

            assert(vf_poles.size() != 0);
            assert(vf_poles.size() == vf_residues.size());

            // compute alpha, lambda first time step; will be 1st order. Note that on first time step only lambda is used (mu and nu are not)
            for (int i = 0; i < vf_poles.size(); i++)           // # poles == # residues
            {
                std::complex<double> x = vf_poles[i];           // pole pk
                std::complex<double> r = vf_residues[i];        // residue rk

                std::complex<double> alpha = std::exp(x * delta);
                std::complex<double> lambda = - (one / x) * ( one + (one - alpha) / (x * delta) );  // 1st order

                res_lambda_sum += r * lambda;
    
                hist_sum += 0.;
                
                alpha_k.push_back(alpha);
                lambda_k.push_back(lambda);
            }

            assert(res_lambda_sum.real() != 0);
            assert(vf_remainder != 0);

            factor = (vf_remainder + 2 * (res_lambda_sum.real()));

            // calculate Gc term
            Gc = (1. - factor) / ( Zref * (1. + factor) );

            // calculate h[n] term; will be 0 on first time step regardless
            h_n = 0.;

            // i_c term is also hence zero on first time step so no need to update from default value

            if (tr_data)
            {
                tr_data << "Time step " << time_step + 1 << '\n';
                tr_data << "\tfactor = " << factor << "\n\tGc = " << Gc << "\n\th_n = " << h_n << "\n\ti_c = " << i_c << '\n'; 
            }
        }

        // calculate constitutive relations for MNA solver

        _y[0].f1 = (1 / Gc);
        _y[0].f0 = (1 / Gc) * i_c;
        _y[0].x = Gc * tr_involts_limited() - i_c;

        assert(converged());
        store_values();
        q_load();   

        _m0.c1 = 1 / _y[0].f1; 
        _m0.c0 = - i_c;
        _m0.x = tr_involts_limited();

        q_accept();             // queue the tr_accept() function
        return converged();
    }

    // here is where we calculate the terms that make up the history term and then propagate to next time step
    void RC_MODEL::tr_accept()
    {
        if (!accepted)
        {
            if (tr_data) tr_data << "Entered tr_accept()...\n";

            dk_n.clear();

            // these SHOULD be the solutions to the current time step
            double vp = tr_involts_limited();               // port voltage from MNA soln (I hope...)
            double ip = Gc * tr_involts_limited() - i_c;    // port current entering block (I hope...)

            // calculate incident wave a[n]
            a_n = ( vp + Zref * ip ) / ( 2 * std::sqrt(Zref) );
            
            // reflected wave b[n] not required, but could be computed here out of interest, if desired!

            // print a[n], a[n-1], and a[n-2] for debugging
            if (tr_data)
            {
                tr_data << "Time step " << time_step + 1 << "\n\t";
                tr_data << "a[" << time_step + 1 << "] = " << a_n << "\ta[" << time_step << "] = " << a_store[0] << "\ta[" << time_step - 1 << "] = " << a_store[1] << '\n'; 
            }

            if (time_step == 0)   // 1st time step
            {
                if (tr_data) tr_data << "\tdk[" << time_step + 1 << "] = \n";

                // calculate current dk[n] term
                for (int i = 0; i < vf_poles.size(); i++)
                {
                    std::complex<double> dk = lambda_k[i] * a_n;    // should be 1st order lambda_k since hasn't been cleared since 1st do_tr()
                    dk_n.push_back(dk);

                    if (tr_data) tr_data << "\tdk[" << time_step + 1 << "](" << i << ") = " << dk << '\n';
                }
            }
            else if (time_step == 1)    // 2nd time step
            {
                if (tr_data) tr_data << "\tdk[" << time_step + 1 << "] = \n";

                // calculate current dk[n] term
                for (int i = 0; i < vf_poles.size(); i++)
                {
                    std::complex<double> dk;
                    dk = alpha_k[i] * dk_store[i] + lambda_k_1[i] * a_n + mu_k_1[i] * a_store[0];
                    dk_n.push_back(dk);

                    if (tr_data) tr_data << "\tdk[" << time_step + 1 << "](" << i << ") = " << dk << '\n';
                }
            }
            else if (time_step > 1) // all other time steps
            {
                if (tr_data) tr_data << "\tdk[" << time_step + 1 << "] = \n";

                // calculate current dk[n] term
                for (int i = 0; i < vf_poles.size(); i++)
                {
                    std::complex<double> dk;
                    dk = alpha_k[i] * dk_store[i] + lambda_k[i] * a_n + mu_k[i] * a_store[0] + nu_k[i] * a_store[1];
                    dk_n.push_back(dk);

                    if (tr_data) tr_data << "\tdk[" << time_step + 1 << "](" << i << ") = " << dk << '\n';
                }
            }

            // store a[n] and dk[n] into previous ones to use for next iteration in tr_advance()
            dk_store = dk_n;                // dk[n]   ----->  dk[n-1]
            a_store[1] = a_store[0];        // a[n-1]  ----->  a[n-2]
            a_store[0] = a_n;               // a[n]    ----->  a[n-1]

            if (tr_data) tr_data << "Exiting tr_accept()...\n";
        }
        accepted = 1;
        if (tr_data) tr_data.close();
    }

    RC_MODEL p1;
    DISPATCHER<CARD>::INSTALL d1(&device_dispatcher, "p|rcm", &p1);
}