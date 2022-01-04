/////Argument Parser for running .C files through root
////Input is of the form:
///       root -l 'file_name.C("-flag1 value -flag2 -flag3 value")'

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <fstream>


TString args;
//Default for run number:
int run;
TString save_as_name = " ";

bool plot_all = false;
bool save = true;
bool no_save = false;
bool plot_cluster_maps = true;
bool plot_strip_info = false;
bool plot_clusters = false;
bool plot_adc = false;
bool nclust_size = false;
int nmax = 100000000;
bool skip_events = false;

void argparser(TString args = "-run 1000"){


//**********************//
//PARSE INPUT ARGUMENTS//
//**********************//
TObjArray *flags = args.Tokenize("-");
int argc = flags->GetEntries();

////////////////////////////////////
//Look for input values after flags//
/////////////////////////////////////
for(int i = 0; i < argc; i++) {
	TString current_val = flags->At(i)->GetName();
	int current_val_length = current_val.Tokenize(" ")->GetEntries();
	TString flag = current_val.Tokenize(" ")->At(0)->GetName();

	if(current_val_length > 1) {
		TString flag_input = current_val.Tokenize(" ")->At(1)->GetName();
		if(flag == "run"){run = flag_input.Atoi();}
		if(flag == "save_as"){save_as_name = flag_input;}
	}
	else{
		if(flag == "plot_all"){plot_all = true;}
		if(flag == "no_save"){no_save = true;}
	}
}

if(plot_all){
    plot_cluster_maps = true;
    plot_strip_info = true;
    plot_clusters = true;
    plot_adc = true;
    nclust_size = true;
  }


}