#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");
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
#include <MAETXIs16.htp>

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
  ntime.allocate("ntime");
  nedades.allocate("nedades");
  minedad.allocate("minedad");
  ntallas.allocate("ntallas");
  mdatos.allocate(1,ntime,1,6,"mdatos");
  Tallas.allocate(1,ntallas,"Tallas");
  Ctot.allocate(1,ntime,1,ntallas,"Ctot");
  msex.allocate(1,ntallas,"msex");
  Wmed.allocate(1,ntallas,"Wmed");
  sigmaR.allocate("sigmaR");
  dt.allocate(1,2,"dt");
  Par_bio.allocate(1,5,"Par_bio");
  hprior.allocate("hprior");
  bprior.allocate("bprior");
 log_Loprior = log(Par_bio(3));
 log_cva_prior = log(Par_bio(4));
 log_b_prior = log(bprior);
  L50prior.allocate("L50prior");
  s1prior.allocate("s1prior");
 log_L50fprior = log(L50prior);
 log_s1prior = log(s1prior);
  nbloques1.allocate("nbloques1");
  ybloques1.allocate(1,nbloques1,"ybloques1");
  nqbloques.allocate("nqbloques");
  yqbloques.allocate(1,nqbloques,"yqbloques");
  opt_qf.allocate("opt_qf");
  opt_bpow.allocate("opt_bpow");
  opt1_fase.allocate("opt1_fase");
  opt_tiposel.allocate("opt_tiposel");
  opt_Lo.allocate("opt_Lo");
  opt_cva.allocate("opt_cva");
  opt_F.allocate("opt_F");
  log_priorRo.allocate("log_priorRo");
  opt_Ro.allocate("opt_Ro");
  opt_devRt.allocate("opt_devRt");
  opt_devNo.allocate("opt_devNo");
  npbr.allocate("npbr");
  pbr.allocate(1,npbr,"pbr");
  ntime_sim.allocate("ntime_sim");
  opt_Frms.allocate("opt_Frms");
  Frms.allocate("Frms");
}

void model_parameters::initializationfunction(void)
{
  log_Lo.set_initial_value(log_Loprior);
  log_cv_edad.set_initial_value(log_cva_prior);
  log_L50.set_initial_value(log_L50fprior);
  log_sigma1.set_initial_value(log_s1prior);
  log_sigma2.set_initial_value(9.2);
  log_b.set_initial_value(log_b_prior);
  log_F.set_initial_value(-1.6);
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
  log_L50.allocate(1,nbloques1,opt1_fase,"log_L50");
  log_sigma1.allocate(1,nbloques1,opt1_fase,"log_sigma1");
  log_sigma2.allocate(1,nbloques1,opt_tiposel,"log_sigma2");
  log_Ro.allocate(opt_Ro,"log_Ro");
  dev_log_Ro.allocate(1,ntime,-10,10,opt_devRt,"dev_log_Ro");
  dev_log_No.allocate(1,nedades,-10,10,opt_devNo,"dev_log_No");
  log_F.allocate(1,ntime,-20,0.7,opt_F,"log_F");
  log_qflo.allocate(1,nqbloques,opt_qf,"log_qflo");
  log_b.allocate(opt_bpow,"log_b");
  log_Lo.allocate(opt_Lo,"log_Lo");
  log_cv_edad.allocate(opt_cva,"log_cv_edad");
  BMflo.allocate(1,ntime,"BMflo");
  #ifndef NO_AD_INITIALIZE
    BMflo.initialize();
  #endif
  Brec.allocate(1,ntime,"Brec");
  #ifndef NO_AD_INITIALIZE
    Brec.initialize();
  #endif
  pred_CPUE.allocate(1,ntime,"pred_CPUE");
  #ifndef NO_AD_INITIALIZE
    pred_CPUE.initialize();
  #endif
  pred_Desemb.allocate(1,ntime,"pred_Desemb");
  #ifndef NO_AD_INITIALIZE
    pred_Desemb.initialize();
  #endif
  likeval.allocate(1,5,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  Neq.allocate(1,nedades,"Neq");
  #ifndef NO_AD_INITIALIZE
    Neq.initialize();
  #endif
  Rpred.allocate(1,ntime,"Rpred");
  #ifndef NO_AD_INITIALIZE
    Rpred.initialize();
  #endif
  Rest.allocate(1,ntime,"Rest");
  Unos_edad.allocate(1,nedades,"Unos_edad");
  #ifndef NO_AD_INITIALIZE
    Unos_edad.initialize();
  #endif
  Unos_anos.allocate(1,ntime,"Unos_anos");
  #ifndef NO_AD_INITIALIZE
    Unos_anos.initialize();
  #endif
  Unos_tallas.allocate(1,ntallas,"Unos_tallas");
  #ifndef NO_AD_INITIALIZE
    Unos_tallas.initialize();
  #endif
  mu_edad.allocate(1,nedades,"mu_edad");
  #ifndef NO_AD_INITIALIZE
    mu_edad.initialize();
  #endif
  sigma_edad.allocate(1,nedades,"sigma_edad");
  #ifndef NO_AD_INITIALIZE
    sigma_edad.initialize();
  #endif
  BDo.allocate(1,ntime,"BDo");
  #ifndef NO_AD_INITIALIZE
    BDo.initialize();
  #endif
  No.allocate(1,nedades,"No");
  #ifndef NO_AD_INITIALIZE
    No.initialize();
  #endif
  prior.allocate(1,7,"prior");
  #ifndef NO_AD_INITIALIZE
    prior.initialize();
  #endif
  yrs.allocate(1,ntime,"yrs");
  #ifndef NO_AD_INITIALIZE
    yrs.initialize();
  #endif
  Desemb.allocate(1,ntime,"Desemb");
  #ifndef NO_AD_INITIALIZE
    Desemb.initialize();
  #endif
  CPUE.allocate(1,ntime,"CPUE");
  #ifndef NO_AD_INITIALIZE
    CPUE.initialize();
  #endif
  cv_index.allocate(1,3,1,ntime,"cv_index");
  #ifndef NO_AD_INITIALIZE
    cv_index.initialize();
  #endif
  nm.allocate(1,ntime,"nm");
  #ifndef NO_AD_INITIALIZE
    nm.initialize();
  #endif
  Lmed_obs.allocate(1,ntime,"Lmed_obs");
  #ifndef NO_AD_INITIALIZE
    Lmed_obs.initialize();
  #endif
  Lmed_pred.allocate(1,ntime,"Lmed_pred");
  #ifndef NO_AD_INITIALIZE
    Lmed_pred.initialize();
  #endif
  edades.allocate(1,nedades,"edades");
  #ifndef NO_AD_INITIALIZE
    edades.initialize();
  #endif
  S1.allocate(1,nbloques1,1,nedades,"S1");
  #ifndef NO_AD_INITIALIZE
    S1.initialize();
  #endif
  S2.allocate(1,nbloques1,1,nedades,"S2");
  #ifndef NO_AD_INITIALIZE
    S2.initialize();
  #endif
  Sel.allocate(1,ntime,1,nedades,"Sel");
  #ifndef NO_AD_INITIALIZE
    Sel.initialize();
  #endif
  F.allocate(1,ntime,1,nedades,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,ntime,1,nedades,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,ntime,1,nedades,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,ntime,1,nedades,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  NM.allocate(1,ntime,1,nedades,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  NMD.allocate(1,ntime,1,ntallas,"NMD");
  #ifndef NO_AD_INITIALIZE
    NMD.initialize();
  #endif
  NDv.allocate(1,ntime,1,ntallas,"NDv");
  #ifndef NO_AD_INITIALIZE
    NDv.initialize();
  #endif
  Nrec.allocate(1,ntime,1,ntallas,"Nrec");
  #ifndef NO_AD_INITIALIZE
    Nrec.initialize();
  #endif
  NVflo.allocate(1,ntime,1,ntallas,"NVflo");
  #ifndef NO_AD_INITIALIZE
    NVflo.initialize();
  #endif
  pred_Ctot.allocate(1,ntime,1,ntallas,"pred_Ctot");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot.initialize();
  #endif
  pred_Ctot_a.allocate(1,ntime,1,nedades,"pred_Ctot_a");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot_a.initialize();
  #endif
  pobs.allocate(1,ntime,1,ntallas,"pobs");
  #ifndef NO_AD_INITIALIZE
    pobs.initialize();
  #endif
  ppred.allocate(1,ntime,1,ntallas,"ppred");
  #ifndef NO_AD_INITIALIZE
    ppred.initialize();
  #endif
  Prob_talla.allocate(1,nedades,1,ntallas,"Prob_talla");
  #ifndef NO_AD_INITIALIZE
    Prob_talla.initialize();
  #endif
  P1.allocate(1,nedades,1,ntallas,"P1");
  #ifndef NO_AD_INITIALIZE
    P1.initialize();
  #endif
  P2.allocate(1,nedades,1,ntallas,"P2");
  #ifndef NO_AD_INITIALIZE
    P2.initialize();
  #endif
  P3.allocate(1,nedades,1,ntallas,"P3");
  #ifndef NO_AD_INITIALIZE
    P3.initialize();
  #endif
  Nv.allocate(1,ntime,1,nedades,"Nv");
  #ifndef NO_AD_INITIALIZE
    Nv.initialize();
  #endif
  NMDv.allocate(1,ntime,1,nedades,"NMDv");
  #ifndef NO_AD_INITIALIZE
    NMDv.initialize();
  #endif
  suma1.allocate("suma1");
  #ifndef NO_AD_INITIALIZE
  suma1.initialize();
  #endif
  suma2.allocate("suma2");
  #ifndef NO_AD_INITIALIZE
  suma2.initialize();
  #endif
  suma3.allocate("suma3");
  #ifndef NO_AD_INITIALIZE
  suma3.initialize();
  #endif
  suma4.allocate("suma4");
  #ifndef NO_AD_INITIALIZE
  suma4.initialize();
  #endif
  penalty.allocate("penalty");
  #ifndef NO_AD_INITIALIZE
  penalty.initialize();
  #endif
  So.allocate("So");
  #ifndef NO_AD_INITIALIZE
  So.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  Linfh.allocate("Linfh");
  #ifndef NO_AD_INITIALIZE
  Linfh.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  h.allocate("h");
  #ifndef NO_AD_INITIALIZE
  h.initialize();
  #endif
  BDp.allocate("BDp");
  #ifndef NO_AD_INITIALIZE
  BDp.initialize();
  #endif
  Npplus.allocate("Npplus");
  #ifndef NO_AD_INITIALIZE
  Npplus.initialize();
  #endif
  Bp_anch.allocate("Bp_anch");
  #ifndef NO_AD_INITIALIZE
  Bp_anch.initialize();
  #endif
  nm1.allocate("nm1");
  #ifndef NO_AD_INITIALIZE
  nm1.initialize();
  #endif
  cuenta1.allocate("cuenta1");
  #ifndef NO_AD_INITIALIZE
  cuenta1.initialize();
  #endif
  Np.allocate(1,nedades,"Np");
  #ifndef NO_AD_INITIALIZE
    Np.initialize();
  #endif
  Zpbr.allocate(1,nedades,"Zpbr");
  #ifndef NO_AD_INITIALIZE
    Zpbr.initialize();
  #endif
  Fpbr.allocate(1,nedades,"Fpbr");
  #ifndef NO_AD_INITIALIZE
    Fpbr.initialize();
  #endif
  Sp.allocate(1,nedades,"Sp");
  #ifndef NO_AD_INITIALIZE
    Sp.initialize();
  #endif
  Bp.allocate(1,npbr,1,ntime_sim,"Bp");
  #ifndef NO_AD_INITIALIZE
    Bp.initialize();
  #endif
  CTPp.allocate(1,nedades,"CTPp");
  #ifndef NO_AD_INITIALIZE
    CTPp.initialize();
  #endif
  Yp.allocate(1,npbr,1,ntime_sim,"Yp");
  #ifndef NO_AD_INITIALIZE
    Yp.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  BD.allocate(1,ntime,"BD");
  BT.allocate(1,ntime,"BT");
  RPR.allocate(1,ntime,"RPR");
  RPRlp.allocate(1,ntime,"RPRlp");
  #ifndef NO_AD_INITIALIZE
    RPRlp.initialize();
  #endif
  SSBo.allocate("SSBo");
  YTPp.allocate(1,npbr,"YTPp");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
   
 yrs=column(mdatos,1);
 Desemb=column(mdatos,2);
 CPUE=column(mdatos,4);
 nm=column(mdatos,6);
 edades.fill_seqadd(minedad,1);
 cv_index(1)=column(mdatos,3);
 cv_index(2)=column(mdatos,5);
 Linf=Par_bio(1);
 k=Par_bio(2);
 M=Par_bio(5);
 
 Unos_edad=1;// lo uso en  operaciones matriciales con la edad
 Unos_anos=1;// lo uso en operaciones matriciales con el a?o
 Unos_tallas=1;// lo uso en operaciones matriciales con el a?o
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{500,2000,5000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-2,1e-5,1e-5}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::userfunction(void)
{
  f =0.0;
 Eval_prob_talla_edad();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_deinteres();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_logverosim();
 Eval_funcion_objetivo();
 if(last_phase()){Eval_CTP();}
}

void model_parameters::Eval_prob_talla_edad(void)
{
 int i, j;
 mu_edad(1)=exp(log_Lo);
 for (i=2;i<=nedades;i++)
  {
  mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
  }
 sigma_edad=exp(log_cv_edad)*mu_edad;
  Prob_talla = ALK( mu_edad, sigma_edad, Tallas);
}

dvar_matrix model_parameters::ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
{
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);
}

void model_parameters::Eval_selectividad(void)
{
 int i,j;
 for (j=1;j<=nbloques1;j++){
 S1(j)=exp(-0.5*square(edades-exp(log_L50(j)))/square(exp(log_sigma1(j))));
    for (i=1;i<=nedades;i++){
      if(edades(i)>=exp(log_L50(j))){
      S1(j,i)= exp(-0.5*square(edades(i)-exp(log_L50(j)))/square(exp(log_sigma2(j))));
      }
 }}
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sel(i)=S1(j);}
       }
   }
}

void model_parameters::Eval_mortalidades(void)
{
 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_edad));
 Z=F+M;
 S=mfexp(-1.0*Z);
}

void model_parameters::Eval_abundancia(void)
{
 int i, j;
 if(opt_Ro<0)
 {
  log_Ro=log_priorRo;
 }
 h=hprior;
 // Biomasa desovante virgen de largo plazo
 No(1)=exp(log_Ro+0.5*square(sigmaR)); //
 for (int j=2;j<=nedades;j++)
 {   No(j)=No(j-1)*exp(-1.*M);}
     No(nedades)=No(nedades)/(1-exp(-1.*M));
 SSBo=sum(elem_prod(No*exp(-dt(1)*M)*Prob_talla,elem_prod(msex,Wmed)));
 alfa=4*h*exp(log_Ro+0.5*square(sigmaR))/(5*h-1);//
 beta=(1-h)*SSBo/(5*h-1);// Reclutamiento
 N(1)=elem_prod(No,exp(dev_log_No));
 BD(1)=sum(elem_prod(elem_prod(N(1),exp(-dt(1)*Z(1)))*Prob_talla,elem_prod(msex,Wmed)));
 Rpred(1)=exp(log_Ro+0.5*square(sigmaR));//
 for (i=1;i<ntime;i++)
 {
    Rpred(i+1)=exp(log_Ro+0.5*square(sigmaR)); // 
    if(i>minedad){
     Rpred(i+1)=(alfa*BD(i-minedad)/(beta+BD(i-minedad)));
    } // 
     N(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i));  // 
     N(i+1)(2,nedades)=++elem_prod(N(i)(1,nedades-1),S(i)(1,nedades-1));
     N(i+1,nedades)+=N(i,nedades)*S(i,nedades);// grupo plus
     BD(i+1)=sum(elem_prod(elem_prod(N(i+1),exp(-dt(1)*Z(i+1)))*Prob_talla,elem_prod(msex,Wmed)));
 }
   Rest=column(N,1);
}

void model_parameters::Eval_deinteres(void)
{
 Nv=N;// solo para empezar los calculos
 for (int i=1;i<ntime;i++)
 {
     Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*M);
     Nv(i+1,nedades)=Nv(i+1,nedades)+Nv(i,nedades)*exp(-1.0*M);// grupo plus
 }
 NDv=elem_prod((Nv*exp(-dt(1)*M))*Prob_talla,outer_prod(Unos_anos,msex));
 BDo=NDv*Wmed;
 RPR=elem_div(BD,BDo);
 RPRlp=BD/SSBo;
}

void model_parameters::Eval_biomasas(void)
{
 NMD=elem_prod(N,mfexp(-dt(1)*Z))*Prob_talla;
 NMD=elem_prod(NMD,outer_prod(Unos_anos,msex));
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel)*Prob_talla;
 BD=NMD*Wmed;
 BMflo=NVflo*Wmed;
 BT=(N*Prob_talla)*Wmed;
}

void model_parameters::Eval_capturas_predichas(void)
{
 pred_Ctot_a=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
 pred_Ctot=pred_Ctot_a*Prob_talla;
 pred_Desemb=pred_Ctot*Wmed;
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));
 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);
}

void model_parameters::Eval_indices(void)
{
   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*pow(BMflo(i),exp(log_b));}
       }
   }
}

void model_parameters::Eval_logverosim(void)
{
 int i;
 suma1=0; suma2=0; penalty=0;
 for (i=1;i<=ntime;i++)
 {
  if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}
 }
}

void model_parameters::Eval_funcion_objetivo(void)
{
 suma3=0; suma4=0; penalty=0;
 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb
 for (int i=1;i<=ntime;i++){
 suma3+=-nm(i)*sum(elem_prod(pobs(i),log(ppred(i))));
 }
 likeval(3)=suma3;//
 if(active(dev_log_Ro)){
 likeval(4)=1./(2*square(sigmaR))*norm2(dev_log_Ro);}
 if(active(dev_log_No)){
 likeval(5)=1./(2*square(sigmaR))*norm2(dev_log_No);}
 if (active(log_F)){
 penalty+=1000*norm2(log_F-mean(log_F));}
 f=(sum(likeval)+penalty);
 if(last_phase){
 f=sum(likeval);}
}

void model_parameters::Eval_CTP(void)
{
  for (int i=1;i<=npbr;i++){ // ciclo de PBR
  Np=N(ntime);
  Sp=S(ntime);
  for (int j=1;j<=ntime_sim;j++){ // ciclo de a?os
  if(j==1){
  Np(1)=(alfa*BD(ntime)/(beta+BD(ntime)));}
  if(j>1){
  Np(1)=(alfa*Bp(i,j-1)/(beta+Bp(i,j-1)));}
  Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
  Np(nedades)+=Np(nedades)*Sp(nedades);
  Fpbr=F(ntime)*pbr(i);//
 if(opt_Frms>0)//agregada 1 Activado Frms/ -1 Activa F ultimo a?o
  {
   Fpbr=Frms*pbr(i);//agregada Activar o desactivar -1
  }
  Zpbr=Fpbr+M;
  Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dt(1)*Zpbr))*Prob_talla,elem_prod(msex,Wmed)));
  CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  Yp(i,j)=sum(elem_prod(CTPp*Prob_talla,Wmed));
  Sp=exp(-1.*Zpbr);
  }
 YTPp(i)=Yp(i,1);//agregada. Ver en STD
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
 report << "years" << endl;
 report << yrs << endl;
 report << "CPUE_obs" << endl;
 report << CPUE << endl;
 report << "CPUE_pred" << endl;
 report << pred_CPUE << endl;
 report << "Desemb_obs" << endl;
 report << Desemb << endl;
 report << "Desemb_pred" << endl;
 report << pred_Desemb << endl;
 report << "Lmed_obs" << endl;
 report << Lmed_obs << endl;
 report << "Lmed_pred" << endl;
 report << Lmed_pred << endl;
 report << "BD" << endl;
 report << BD << endl;
 report << "BT" << endl;
 report << BT << endl;
 report << "BV" << endl;
 report << BMflo << endl;
 report << "R_pred" << endl;
 report << Rpred<< endl;
 report << "R_Est" << endl;
 report << column(N,1)<< endl;
 report << "F " << endl;
 report << exp(log_F) << endl;
 report << "Edades"<< endl;
 report << edades<< endl;
 report <<"N"<<endl;
 report <<N<< endl;
 report <<"Sel_f"<<endl;
 report <<Sel<< endl;
 report <<"pobs"<< endl;
 report <<pobs<< endl;
 report <<"ppred"<< endl;
 report <<ppred<< endl;
 report << "Tallas"<< endl;
 report << Tallas<< endl;
 report << "Prob_talla" << endl;
 report << Prob_talla << endl;
 report << "BDo" << endl;
 report << SSBo << endl;
 report << "Lmed" << endl;
 report << mu_edad<< endl;
 report << "likeval"<<endl;
 report << likeval << endl;
 report << "Ro"<<endl;
 report << log_Ro << endl;
}

void model_parameters::final_calcs()
{
 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

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
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 
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
