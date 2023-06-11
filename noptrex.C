// includes and paths
#include "TCanvas.h"
#include "TLeaf.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TTreeFormula.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include <sys/stat.h>

#define CONF "./meta/run%05d.conf"
#define RUN  "./data/run%05d.bin"
#define SPIN "./spin/run%05d.spin"
#define MAXCH 256
#define LINE_LENGTH 4096
#define MAX_LEAFS 1024
#define DT 5000000
#define ET 10000


TTree* noptrex(int run=-1, int verbose=0, TString opt="")
{

  
  // parse config file to environment variables with same name
  #define MAXB 4096 // max number of characters per input line, max length of option array
    char charbuf[MAXB], *strptr; // general purpose strings for parsing
    FILE* fp=fopen(Form(CONF,run),"r");  if (fp) while (fgets(charbuf,MAXB,fp)) {
        if ((strptr=strchr(charbuf,'#'))) *strptr='\0'; // get rid of comments
        if ((strptr=strchr(charbuf,'='))) *strptr='\0', setenv(charbuf,strptr+1,1);
      } else {
      printf("Config file %s not found\n", Form(CONF,run));  return 0;
    } fclose(fp);

  // macros for parsing option variables from environment
  #define INTOPT(NAME,DEF) NAME = getenv((#NAME)) ? atoi(getenv((#NAME))) : DEF
  #define STROPT(NAME,DEF) NAME[MAXB]=""; strcpy(NAME, getenv(#NAME) ? getenv(#NAME) : DEF)
  #define ARROPT(NAME,DEF) NAME[MAXB]={0}; int NAME##_n=0;		\
  strcpy(charbuf, getenv(#NAME)?getenv(#NAME):DEF ); strptr=strtok(charbuf," "); \
  while (strptr) NAME[NAME##_n++]=atoi(strptr), strptr=strtok(0," ");

  // option variables, set from the environment by above macros
  int INTOPT(decimation, 6);	// decimation setting in hardware
  int INTOPT(raw, 0);		// 1:record the raw waveforms before the windows, 0: record only windows
  int INTOPT(nwin,  2 );	// number of decimation windows
  int ARROPT(wchn, "1 1");	// channel number  of each window
  int ARROPT(wbeg, "100 1000"); // first sample    of each window
  int ARROPT(wend, "200 2000"); // last sample+1   of each window
  int ARROPT(wsum, "10 100");   // #samples to add in each window
  int ARROPT(wdec, "2 5");      // #bits decimated in each window
  int INTOPT(npulses,6000);     // #pulses in data file
  uint32_t INTOPT(mask, 15);    // LSB (rightmost bit) is channel zero, followed by channel 1, 2..7, read right to left.
  uint32_t INTOPT(reclength, 4096);
  uint32_t INTOPT(posttrig, 80);
  uint32_t ARROPT(DCoffset, "32768 32768 32768 32768 32768 32768 32768 32768");


  int ARROPT(chsf,"4 5");	// window numbers
  int ARROPT(tlosf,"0 0");	// range of pickup loops to integrate
  int ARROPT(thisf,"750 750");
  int ARROPT(vlosf,"8590 8420 8525 8700 8540 7900 8100 8240 8590 8700 8505 8200 8900 8375 8505 7900 8540 8375"); // vlosf[nseq+1][chsf_n] range of thresholds for each spin state in groups of (chsf), last group is continuation thresholds (no pickup loop signal)
  int ARROPT(vhisf,"8625 8520 8540 8840 8590 8100 8300 8520 8625 8840 8525 8375 9100 8420 8525 8100 8590 8420");


//////////////////////////////////////////////////////////////////////////////////////////



  // TFile * hFile = hFile = TFile::Open(Form("run%d_raw.root",run),"recreate");  //create new root file for run
  FILE * fFile = fFile = fopen(Form(RUN,run),"r");  //open data file for run
  FILE * sFile = sFile = fopen(Form(SPIN,run),"r");  //open data file for run

  //calculate buffer size for each window and add to array
  int w_size[nwin];
  for (int i = 0; i<nwin; i++){
    w_size[i]=(wend[i]-wbeg[i])/wsum[i];
  }

  //create a buffer for each window
  UInt_t h[5];
  UShort_t w0[w_size[0]];
  UShort_t w1[w_size[1]];
  UShort_t w2[w_size[2]];
  UShort_t w3[w_size[3]];
  UShort_t w4[w_size[4]];
  UShort_t w5[w_size[5]];
  UShort_t w6[w_size[6]];
  UShort_t w7[w_size[7]];
  UShort_t w8[w_size[8]];

  UShort_t * w [9] = {w0,w1,w2,w3,w4,w5,w6,w7,w8};  //create array containing pointers to window buffers

  TTree* t = new TTree( Form("run%d_raw",run),Form("raw data for run %d",run) );


  //create branch for each window
  t->Branch("h",&h,"h[5]/i");
  for (int i=0;i<nwin;i++){
    t->Branch(Form("w%d",i),w[i],Form("w%d[%d]/s",i,w_size[i]));
  }

  //read data into buffer for each window and fill the TTree
  for (int ev=0; ev < npulses; ++ev) {
    fread(h,sizeof(h),1,fFile);
    for (int i=0;i<nwin;i++){
      fread(w[i],w_size[i]*2,1,fFile);
    }


    t->Fill();
  }

/////////////////////////////////////////////////////////

int state[]={-1,-1,-1,-1};
unsigned int clock=0;
int states[100000][4];		// buffered output to revise before writing





t->Branch("pulse",&state[0],"pulse/I");
t->Branch("step",&state[1],"step/I");
t->Branch("sequence",&state[2],"sequence/I");
t->Branch("resets",&state[3],"resets/I");


TTreeFormula *th=new TTreeFormula("","h",t), *tf[chsf_n];
  for (int ch=0; ch<chsf_n; ++ch) {
    tf[ch]=new TTreeFormula("",Form("w%d",chsf[ch]),t);
  }
  long int* signal=new long int[chsf_n];

  for (int ev=0; ev < npulses; ++ev) {
    
    // average signals over specified ranges
    t->GetEntry(ev); 
    unsigned int tclock=th->EvalInstance(3);
    if (abs(int(tclock-clock)-DT)>ET && abs(int((1<<31)-(clock-tclock))-DT)>ET) {
      if (state[2]>=0) for (int rev=ev; states[--rev][2]==state[2]; ) states[rev][2]=-1;
      state[0]=-1; state[1]=-1; state[2]=-1; state[3]++;
      if (state[3]) printf("skip ");
    }
    clock = tclock;
    for (int ch=0; ch<chsf_n; ++ch) {
      signal[ch]=0;
      for (int t=tlosf[ch]; t<thisf[ch]; ++t) {
	    signal[ch]+=tf[ch]->EvalInstance(t);
      }
      signal[ch] *= ( (1<<wdec[chsf[ch]]) / float(thisf[ch]-tlosf[ch]) / wsum[chsf[ch]] );
    }
    
    // match average signals to thresholds
    int seq=vlosf_n/chsf_n;
    while (--seq>=0) {
      bool check=true;
      for (int ch=0; ch<chsf_n; ++ch) {
      	if ( signal[ch]<vlosf[seq*2+ch] ||
	        signal[ch]>vhisf[seq*2+ch] ) check=false;
      }
      if (check) break;
    }
    if (seq==vlosf_n/chsf_n-1) state[0]++;
    else state[0]=0, state[1]=seq, state[2]+=!seq;

    for (int i=0; i<4; ++i) states[ev][i]=state[i];


    t->GetBranch("pulse")->Fill();
    t->GetBranch("step")->Fill();
    t->GetBranch("sequence")->Fill();
    t->GetBranch("resets")->Fill();
  }




  fclose(fFile);
  fclose(sFile);
  // t->Write("",TObject::kOverwrite);
  // t1->Write("",TObject::kOverwrite);
  // t->Write();
  // t1->Write();
  // hFile->Write();
  // hFile->Close();

  // c.Stop();
  // c.Print();


  //////////////////////////////////////////////////////////////////////////////////////////

  return t;
}
