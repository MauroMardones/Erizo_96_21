#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <PBRmodbento.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  M.allocate("M");
  h.allocate("h");
  dt.allocate("dt");
  npbr.allocate("npbr");
  ratio.allocate(1,npbr,"ratio");
  nedades.allocate("nedades");
  Msex.allocate(1,nedades,"Msex");
  Wm.allocate(1,nedades,"Wm");
  nyrs.allocate("nyrs");
  Sel.allocate(1,nyrs,1,nedades,"Sel");
}

void model_parameters::initializationfunction(void)
{
  log_Frms.set_initial_value(-0.69);
  log_Fspr.set_initial_value(-0.69);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_Frms.allocate(1,nyrs,"log_Frms");
  log_Fspr.allocate(1,npbr,1,nyrs,"log_Fspr");
  N.allocate(1,nedades,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  Fcr.allocate(1,500,"Fcr");
  #ifndef NO_AD_INITIALIZE
    Fcr.initialize();
  #endif
  Z.allocate(1,nedades,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(1,nedades,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(1,nedades,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  C.allocate(1,nedades,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  BD.allocate(1,500,"BD");
  #ifndef NO_AD_INITIALIZE
    BD.initialize();
  #endif
  Y.allocate(1,500,"Y");
  #ifndef NO_AD_INITIALIZE
    Y.initialize();
  #endif
  BDLP.allocate(1,500,"BDLP");
  #ifndef NO_AD_INITIALIZE
    BDLP.initialize();
  #endif
  YLP.allocate(1,500,"YLP");
  #ifndef NO_AD_INITIALIZE
    YLP.initialize();
  #endif
  RLP.allocate(1,500,"RLP");
  #ifndef NO_AD_INITIALIZE
    RLP.initialize();
  #endif
  ratio_pbr.allocate(1,npbr,"ratio_pbr");
  #ifndef NO_AD_INITIALIZE
    ratio_pbr.initialize();
  #endif
  aux.allocate("aux");
  #ifndef NO_AD_INITIALIZE
  aux.initialize();
  #endif
  MRS.allocate("MRS");
  #ifndef NO_AD_INITIALIZE
  MRS.initialize();
  #endif
  Bmrs.allocate("Bmrs");
  #ifndef NO_AD_INITIALIZE
  Bmrs.initialize();
  #endif
  Fmrs.allocate("Fmrs");
  #ifndef NO_AD_INITIALIZE
  Fmrs.initialize();
  #endif
  Bo.allocate("Bo");
  #ifndef NO_AD_INITIALIZE
  Bo.initialize();
  #endif
  No.allocate(1,nedades,"No");
  #ifndef NO_AD_INITIALIZE
    No.initialize();
  #endif
  BDo.allocate("BDo");
  #ifndef NO_AD_INITIALIZE
  BDo.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  BPR.allocate("BPR");
  #ifndef NO_AD_INITIALIZE
  BPR.initialize();
  #endif
  YPR.allocate("YPR");
  #ifndef NO_AD_INITIALIZE
  YPR.initialize();
  #endif
  BPRLP.allocate(1,nyrs,"BPRLP");
  #ifndef NO_AD_INITIALIZE
    BPRLP.initialize();
  #endif
  Yopt.allocate(1,nyrs,"Yopt");
  #ifndef NO_AD_INITIALIZE
    Yopt.initialize();
  #endif
  BPR2.allocate("BPR2");
  #ifndef NO_AD_INITIALIZE
  BPR2.initialize();
  #endif
  BPRLP2.allocate("BPRLP2");
  #ifndef NO_AD_INITIALIZE
  BPRLP2.initialize();
  #endif
  penalty.allocate("penalty");
  #ifndef NO_AD_INITIALIZE
  penalty.initialize();
  #endif
  dy01.allocate("dy01");
  #ifndef NO_AD_INITIALIZE
  dy01.initialize();
  #endif
  F01.allocate("F01");
  #ifndef NO_AD_INITIALIZE
  F01.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
 for (int k=1;k<=nyrs;k++){
  N(1)=1.;
  for (int j=2;j<=nedades;j++)
  { N(j)=N(j-1)*exp(-1.*M);
    N(nedades)=N(nedades)/(1-exp(-1.*M));
  }
  Bo=sum(elem_prod(elem_prod(N*exp(-dt*M),Msex),Wm));
  alfa=4*h/(5*h-1);
  beta=(1-h)/(5*h-1)*Bo;
 //-estima MRS---------------------------------------------
  F=exp(log_Frms(k))*Sel(k);
  Z=F+M;
  S=exp(-1.*Z);
  // se estima la sobrevivencia por edad y a?o
  for (int i=2;i<=nedades;i++){
  N(i)=N(i-1)*exp(-Z(i-1));
  }
  N(nedades)=N(nedades)/(1-exp(-Z(nedades)));
  BPR=sum(elem_prod(elem_prod(elem_prod(N,exp(-dt*Z)),Msex),Wm));
  YPR=sum(elem_prod(elem_prod(elem_div(F,Z),elem_prod(1.-S,N)),Wm));
  BPRLP(k)=alfa*BPR-beta;//Biomasa de equilibrio
  Yopt(k)=YPR*(alfa*BPRLP(k)/(beta+BPRLP(k)));// rendimiento de equilibrio
 //--Estima BRP ratio y F%---------------------------------
  for (int j=1;j<=npbr;j++){
    Z=exp(log_Fspr(j,k))*Sel(k)+M;
    for (int i=2;i<=nedades;i++){
     N(i)=N(i-1)*exp(-Z(i-1));
    }
  N(nedades)=N(nedades)/(1-exp(-Z(nedades)));
  BPR2=sum(elem_prod(elem_prod(elem_prod(N,exp(-dt*Z)),Msex),Wm));
  BPRLP2=alfa*BPR2-beta;// Biomasa de equilibrio
  ratio_pbr(j)=BPRLP2/Bo;}
  penalty=10000*norm2(ratio_pbr-ratio);
  f+=-Yopt(k)+penalty;
  }
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
   //
  Fcr.fill_seqadd(0,0.02);
  for (int j=1;j<=500;j++) // ciclo de F's
  {
  F=Fcr(j)*Sel(nyrs);
  Z=F+M;
  S=exp(-1.*Z);
  N(1)=1;
  // se estima la sobrevivencia por edad y a?o
     for (int i=2;i<=nedades;i++){ // ciclo de a?os
     N(i)=N(i-1)*exp(-Z(i-1));
     }
     N(nedades)=N(nedades)/(1-exp(-Z(nedades)));
     C=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
     Y(j)=sum(elem_prod(C,Wm));// rendimiento 
     BD(j)=sum(elem_prod(elem_prod(elem_prod(N,exp(-dt*Z)),Msex),Wm));
  }
  dy01=0.1*(Y(2)-Y(1))/Fcr(2);
  for (int j=1;j<=500-1;j++) // ciclo de F's
  {
    if((Y(j+1)-Y(j))/Fcr(2)-dy01<0){
    F01=0.5*(Fcr(j)+Fcr(j+1));j=500;}
   }   
     BDLP=alfa*BD-beta;//
     RLP=elem_div(alfa*BDLP,beta+BDLP);
     YLP=elem_prod(Y,RLP);
   int i;
   i=1;
   aux=1;
  report << "Biological Reference Points -BRPmodel- " << endl;
  report << "IFOP 2014 (ccr) v2.0 " << endl;
  report << "h" << endl;
  report << h << endl;
  report << "SSBPRo " << endl;
  report << Bo << endl;
  report << "Fmsy" << endl;
  report <<  exp(log_Frms) << endl;
  report << "MSY/R" << endl;
  report <<  Yopt << endl;
  report << "SSBmsy/R" << endl;
  report <<  BPRLP << endl;
  report << "SSBmsy/SSBo" << endl;
  report << BPRLP/Bo<< endl;
  report << "%SSBo"<<endl;
  report << ratio_pbr<<endl; 
  report << "F(%SSBo)"<<endl;
  report << exp(log_Fspr)<<endl; 
  report << "Fo1"<<endl;
  report << F01 <<endl;
  report << "Fcr_YPR_SSBPR_SSBeq_Req_SSBeq/SSBo_Yeq" << endl;
   while(aux==1){
   report << Fcr(i) <<" "<<Y(i)<<" "<<BD(i)<<" "<<BDLP(i) <<" "<<RLP(i)<<"  "<<BDLP(i)/BDLP(1)<<" " <<YLP(i)<<endl;
   aux=1;
   i=i+1;
   if(YLP(i)<0.001*max(YLP)){
   exit(1);} 
    }
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
