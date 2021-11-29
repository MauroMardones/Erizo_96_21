DATA_SECTION

  init_number M  
  init_number h
  init_number dt
  init_int npbr
  init_vector ratio(1,npbr)
  init_int nedades
  init_vector Msex(1,nedades)
  init_vector Wm(1,nedades)
  init_int nyrs
  init_matrix Sel(1,nyrs,1,nedades)

INITIALIZATION_SECTION

  log_Frms        -0.69
  log_Fspr        -0.69

PARAMETER_SECTION

 init_vector log_Frms(1,nyrs) 
 init_matrix log_Fspr(1,npbr,1,nyrs)


//Defino las variables de estado 
 vector N(1,nedades)
 vector Fcr(1,500)
 vector Z(1,nedades)
 vector F(1,nedades)
 vector S(1,nedades)
 vector C(1,nedades)
 vector BD(1,500)
 vector Y(1,500)
 vector BDLP(1,500)
 vector YLP(1,500)
 vector RLP(1,500)
 vector ratio_pbr(1,npbr)

 number aux
 number MRS
 number Bmrs
 number Fmrs
 number Bo
 vector No(1,nedades)


 number BDo
 number alfa
 number beta
 number BPR
 number YPR
 vector BPRLP(1,nyrs)
 vector Yopt(1,nyrs)

 number BPR2
 number BPRLP2
 number penalty
 number dy01
 number F01


 objective_function_value f


PROCEDURE_SECTION


 for (int k=1;k<=nyrs;k++){

// condici?n inicial Bo
  N(1)=1.;
  for (int j=2;j<=nedades;j++)
  { N(j)=N(j-1)*exp(-1.*M);
    N(nedades)=N(nedades)/(1-exp(-1.*M));
  }
  Bo=sum(elem_prod(elem_prod(N*exp(-dt*M),Msex),Wm));

// parametros S/R
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




REPORT_SECTION
   //


//vector de mortalidad por pesca a evaluar
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
    
