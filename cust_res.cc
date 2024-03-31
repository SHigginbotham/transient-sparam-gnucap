/*
  SIMPLE TEST DEVICE MODEL PLUGIN FOR GNUCAP
  IMPLEMENTS V = I * R*R

  Created by:     Sean Higginbotham for M.A.I project
                  Supervisor: Dr. Justin King
                  Department of Electronic and Electrical Engineering,
                  Trinity College Dublin, Ireland, September 2023 - April 2024
*/

#include "globals.h"
#include "e_elemnt.h"

namespace
{
  class CUST_RESISTANCE : public ELEMENT
  {
  private:
		// copy constructor used to create instances of the device from the static instance
		explicit CUST_RESISTANCE(const CUST_RESISTANCE &p) : ELEMENT(p) {}
		// should use to de-allocate any dynamic memory if create any; otherwise not needed
		// ~CUST_RESISTANCE() {}
  public:
		// default constructor used by dispatcher to create static instance
		explicit CUST_RESISTANCE() : ELEMENT() {}

  private:
		// parser functions
		char id_letter() const { return 'a'; }              // id used by parser
		std::string value_name() const { return "r"; }      // value 'name'
		std::string dev_type() const { return "cust_res"; } // model 'type'
		int max_nodes() const { return 2; }                 // max # of ports parser will allow
		int min_nodes() const { return 2; }                 // min # of ports parser will allow
		// allocation functions
		int net_nodes() const { return 2; }    // actual # of ports for MNA matrix
		int matrix_nodes() const { return 2; } // the # of nodes that will be stamped into matrix
		// port functions
		std::string port_name(int i) const
		{
				// must be 0 (+ve) or 1 (-ve) terminal
				assert(i >= 0);
				assert(i < 2);
				static std::string names[] = {"p", "n"};
				return names[i];
		}

		// used to copy this copy this instance
		CARD *clone() const { return new CUST_RESISTANCE(*this); }

		// for ELEMENT derived classes, need to override the following functions
		// relating to TR and AC analysis:
		// ----------------------------------------------------------------------
		// notifies the admittance matrix and LU matrix of the nodes used by this device.
		// Call *_passive(), *_active(), or *_extended() depending on device type
		void tr_iwant_matrix() { tr_iwant_matrix_passive(); }
		// obtain port voltages
		double tr_involts() const { return tr_outvolts(); }
		double tr_involts_limited() const { return tr_outvolts_limited(); }
		// same as TR version but for AC
		COMPLEX ac_involts() const { return ac_outvolts(); }
		void ac_iwant_matrix() { tr_iwant_matrix_passive(); }

		// transient analysis functions
		void tr_load() { tr_load_passive(); }     // load amittance matrix and current vector with values calculated during do_tr
		void tr_unload() { tr_unload_passive(); } // removes component from MNA matrix
		void tr_begin();
		bool do_tr();
		// there are many other functions (for TR analysis) but they are evidently
		// not required to be overriden unless explicitly required.
  };

  // use the following notation
  /*
		use y.x = i, y.f0 = 0, y.f1 = ohms^2
		use m.x = v, m.c0 = 0, m.c1 = mhos^2
  */

  // is called at beginning of analysis; sets up params and such.
  // main purpose is to just initialise the FPOLY and CPOLY
  void CUST_RESISTANCE::tr_begin()
  {
		// need to call base class virtual first
		ELEMENT::tr_begin();

		// values that will remain constant during sim
		_y[0].f0 = _m0.c0 = 0;
		_y[0].f1 = (value() != 0.) ? value() * value() : OPT::shortckt; // R^2
		_m0.c1 = 1. / _y[0].f1;                                         // 1 / R^2
  }

  // does most of the 'real work'. Optionally calls tr_eval() to do the stuff
  // Can forego and do own stuff too!
  bool CUST_RESISTANCE::do_tr()
  {
		// (1) assign non-constant FPOLY terms
		_y[0].x = tr_involts_limited() * (1. / _y[0].f1); // i[n]

		// (3) necessary process functions
		set_converged(conv_check()); // check convergence
		store_values();              // push values to previous iteration (i.e., store them!)
		q_load();                    // queue for adding to MNA matrix

		// (3) 'convert' to CPOLY values
		_m0.x = tr_involts_limited(); // v[n]

		return converged();
  }

  CUST_RESISTANCE p1;
  DISPATCHER<CARD>::INSTALL d2(&device_dispatcher, "p|custr", &p1);
}