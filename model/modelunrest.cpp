/*
 * modelunrest.cpp
 *
 *  Created on: 24/05/2016
 *      Author: mdw2
 */

#include "modelunrest.h"
#include <stdlib.h>
#include <string.h>

ModelUnrest::ModelUnrest(PhyloTree *tree, string model_params)
	: ModelMarkov(tree, false)
{
	num_params = getNumRateEntries() - 1;
	model_parameters = new double [num_params];
	for (int i=0; i< num_params; i++) model_parameters[i] = 1;
	setRates();
	if (model_params != "") {
		//cout << "WARNING: Supplying model params to constructor not yet properly implemented -- ignored" << endl;
		// TODO: parse model_params into model_parameters, then call setRates().
		int end_pos = 0;
		cout << __func__ << " " << model_params << endl;
		for (int i = 0; i < 12; i++) {
			int new_end_pos;
			try {
				rates[i] = convert_double(model_params.substr(end_pos).c_str(), new_end_pos);
			} catch (string &model_params) {
				outError(model_params);
			}
		
			end_pos += new_end_pos;
			if (rates[i] <= 0.0)
				outError("Non-positive rates found");
			if (i == 11 && end_pos < model_params.length())
				outError("String too long ", model_params);
			if (i < 11 && end_pos >= model_params.length())
				outError("Unexpected end of string ", model_params);
			if (end_pos < model_params.length() && model_params[end_pos] != ',')
				outError("Comma to separate rates not found in ", model_params);
			end_pos++;
		}
		num_params = 0;
		writeInfo(cout);
	}
    name = "UNREST";
    full_name = "Unrestricted model (non-reversible)";
    ModelMarkov::init(FREQ_ESTIMATE);
}

/* static */ bool ModelUnrest::validModelName(string model_name) {
	return (model_name == "UNREST");
}

void ModelUnrest::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int i, ndim = getNDim();

	for (i = 1; i <= ndim; i++) {
		lower_bound[i] = 0.01;
		upper_bound[i] = 100.0;
		bound_check[i] = false;
	}
}

/*
 * Set rates from model_parameters
 */
void ModelUnrest::setRates() {
	// For UNREST, parameters are simply the off-diagonal rate matrix entries
	// (except [4,3] = rates[11], which is constrained to be 1)
	memcpy(rates, model_parameters, num_params*sizeof(double));
	rates[num_params]=1;
}

// void ModelGTR::readParameters(const char *file_name) { 
// 	try {
// 		ifstream in(file_name);
// 		if (in.fail()) {
// 			outError("Invalid model name ", file_name);
//         }
// 		cout << "Reading model parameters from file " << file_name << endl;
// 		readRates(in);
// 		readStateFreq(in);
// 		in.close();
// 	}
// 	catch (const char *str) {
// 		outError(str);
// 	} 
// 	num_params = 0;
// 	writeInfo(cout);
// }
// 
// void ModelGTR::readRates(string str) throw(const char*) {
// 	int nrates = getNumRateEntries();
// 	int end_pos = 0;
// 	cout << __func__ << " " << str << endl;
// 	if (str.find("equalrate") != string::npos) {
// 		for (int i = 0; i < nrates; i++)
// 			rates[i] = 1.0;
// 	} else for (int i = 0; i < nrates; i++) {
// 		int new_end_pos;
// 		try {
// 			rates[i] = convert_double(str.substr(end_pos).c_str(), new_end_pos);
// 		} catch (string &str) {
// 			outError(str);
// 		}
// 		end_pos += new_end_pos;
// 		if (rates[i] <= 0.0)
// 			outError("Non-positive rates found");
// 		if (i == nrates-1 && end_pos < str.length())
// 			outError("String too long ", str);
// 		if (i < nrates-1 && end_pos >= str.length())
// 			outError("Unexpected end of string ", str);
// 		if (end_pos < str.length() && str[end_pos] != ',')
// 			outError("Comma to separate rates not found in ", str);
// 		end_pos++;
// 	}
// 	num_params = 0;
// }

void ModelUnrest::writeInfo(ostream &out) {
	if (num_states == 4) {
		out << "Rate parameters:";
		//out.precision(3);
		//out << fixed;
		out << "  A-C: " << rates[0];
		out << "  A-G: " << rates[1];
		out << "  A-T: " << rates[2];
		out << "  C-A: " << rates[3];
		out << "  C-G: " << rates[4];
		out << "  C-T: " << rates[5];
		out << "  G-A: " << rates[6];
		out << "  G-C: " << rates[7];
		out << "  G-T: " << rates[8];
		out << "  T-A: " << rates[9];
		out << "  T-C: " << rates[10];
		out << "  T-G: " << rates[11];
		out << endl;
	}
}