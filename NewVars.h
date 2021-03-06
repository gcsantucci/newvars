//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 28 15:10:55 2016 by ROOT version 6.06/02
// from TChain h1/
//////////////////////////////////////////////////////////

#ifndef NewVars_h
#define NewVars_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class NewVars {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nring;
   UInt_t          nrun;
   Int_t           nev;
   Int_t           nsub;
   UChar_t         cate;
   Float_t         potot;
   UShort_t        nhit;
   Float_t         pomax;
   Float_t         potota;
   UShort_t        nhita;
   UShort_t        nhitac;
   Float_t         pomaxa;
   Float_t         wall;
   Float_t         evis;
   Float_t         rtsum;
   Float_t         rtmax;
   Float_t         wlen;
   UChar_t         ip[10];   //[nring]
   Float_t         pos[3];
   Float_t         dir[10][3];   //[nring]
   Float_t         dirtot[3];
   Float_t         ang[10];   //[nring]
   Float_t         rtot[10];   //[nring]
   Float_t         amom[10];   //[nring]
   Float_t         rtote[10];   //[nring]
   Float_t         amome[10];   //[nring]
   Float_t         rtotm[10];   //[nring]
   Float_t         amomm[10];   //[nring]
   UChar_t         nsube;
   UChar_t         ndcy;
   UChar_t         ngate;
   UChar_t         nbye;
   Float_t         probms[10][6];   //[nring]
   Float_t         prmslg[10][6];   //[nring]
   Int_t           date[3];
   Int_t           time[4];
   Float_t         elpsday;
   Int_t           numpo[10];   //[nring]
   Float_t         apos[3];
   Float_t         adir[3];
   Float_t         aang;
   Float_t         agood;
   Float_t         wgain;
   Int_t           nbad;
   Int_t           nbada;
   Float_t         msdir[10][3][6];   //[nring]
   Float_t         amomp[10];
   Float_t         ange[10];
   Float_t         angm[10];
   Float_t         angp[10];
   Int_t           ntot[10];
   Float_t         probth[10][6];
   Float_t         probpt[10][6];
   Float_t         pi0like[2];
   Float_t         pi0_e[2][2];
   Float_t         pi0_dir[2][2][3];
   Float_t         pi0mass[2];
   Float_t         evisold;
   Float_t         evisoldxe;
   Float_t         evisnew;
   Int_t           nmue;
   UChar_t         etype[10];   //[nmue]
   Float_t         etime[10];   //[nmue]
   Float_t         epos[10][3];   //[nmue]
   Float_t         edir[10][3];   //[nmue]
   Float_t         egood[10];   //[nmue]
   Float_t         ehit[10];   //[nmue]
   Float_t         mueprob[2];
   Int_t           nmue_sel;
   Int_t           Rnring;
   Float_t         Rdir[30][3];   //[Rnring]
   Float_t         Rang[30];   //[Rnring]
   UChar_t         Riring;
   Float_t         Rtwout[30];   //[Rnring]
   Float_t         Rtwith[30];   //[Rnring]
   Float_t         Alwout;
   Float_t         Alwith;
   Float_t         Qsmi;
   Float_t         Qsmo;
   Float_t         Qexi;
   Float_t         Qexo;
   Float_t         Pe5d;
   Float_t         En5d;
   Float_t         Eh5d;
   Float_t         Pe5do;
   Float_t         En5do;
   Float_t         Eh5do;
   Float_t         Rtadd;
   Float_t         Pdgeta;
   Float_t         Pd5d;
   Float_t         Pdthre;
   Float_t         Pd5do;
   Float_t         Delpd;
   Float_t         Ropena[30];   //[Rnring]
   Int_t           Maxth;
   Float_t         Pkang;
   Float_t         Qrfct;
   Float_t         Pdfct;
   Float_t         Pkfct;
   Float_t         Agfct;
   Float_t         Dlfct;
   Int_t           Iflag;
   Float_t         Pmfct;
   Float_t         Imfct;
   Float_t         Rilike;
   Int_t           ri_ver;
   Float_t         ri_pid;
   Int_t           ri_nring;
   Float_t         ri_flag[10];   //[ri_nring]
   Float_t         ri_dlfct[10];   //[ri_nring]
   Float_t         ri_pdfct[10];   //[ri_nring]
   Float_t         ri_pkfct[10];   //[ri_nring]
   Float_t         ri_vafct[10];   //[ri_nring]
   Float_t         ri_total[10];   //[ri_nring]
   Float_t         ri_dir[10][3];   //[ri_nring]
   Float_t         ri_imfct[10];   //[ri_nring]
   Float_t         ri_pmfct[10];   //[ri_nring]
   Int_t           npar;
   Float_t         wallv;
   UChar_t         ipv[50];   //[npar]
   Float_t         posv[3];
   Float_t         dirv[50][3];   //[npar]
   Float_t         pmomv[50];   //[npar]
   Int_t           light_flag[50];   //[npar]
   Int_t           npar2;
   Float_t         wallv2[50];   //[npar2]
   UChar_t         ipv2[50];   //[npar2]
   UChar_t         iorg[50];   //[npar2]
   Float_t         posv2[50][3];   //[npar2]
   Float_t         dirv2[50][3];   //[npar2]
   Float_t         pmomv2[50];   //[npar2]
   Int_t           numnu;
   Char_t          mode;
   Int_t           ipnu[50];   //[numnu]
   Float_t         pnu[50];   //[numnu]
   Float_t         dirnu[50][3];   //[numnu]
   Float_t         flxg[3];
   Float_t         flxh01[3];
   Float_t         kflux[4];
   Float_t         bs71[3];
   Float_t         bs74[3];
   Float_t         flxf[3];
   Float_t         flxh1d[3];
   Float_t         flxb03[3];
   Float_t         flxf03[3];
   Float_t         flxh06[3];
   Float_t         flxh11[3];
   Int_t           scan[2];
   Int_t           nuceff_version;
   Int_t           charge_exchange;
   Int_t           absorbed;
   Int_t           multipi_gen;
   Int_t           scattering;
   Int_t           piless_dcy;
   Int_t           nscndprt;
   Int_t           itrkscnd[1000];   //[nscndprt]
   Int_t           istakscnd[1000];   //[nscndprt]
   Float_t         vtxscnd[1000][3];   //[nscndprt]
   Float_t         pscnd[1000][3];   //[nscndprt]
   Int_t           iprtscnd[1000];   //[nscndprt]
   Float_t         tscnd[1000];   //[nscndprt]
   Int_t           iprntprt[1000];   //[nscndprt]
   Int_t           lmecscnd[1000];   //[nscndprt]
   Int_t           iprnttrk[1000];   //[nscndprt]
   Int_t           iorgprt[1000];   //[nscndprt]
   Float_t         pprnt[1000][3];   //[nscndprt]
   Int_t           iflgscnd[1000];   //[nscndprt]
   Float_t         pprntinit[1000][3];   //[nscndprt]
   Float_t         vtxprnt[1000][3];   //[nscndprt]
   Int_t           iprntidx[1000];   //[nscndprt]
   Int_t           nchilds[1000];   //[nscndprt]
   Int_t           ichildidx[1000];   //[nscndprt]
   Int_t           Npvcscnd;
   Int_t           Nchildsvc[100];   //[Npvcscnd]
   Int_t           Ichildidxvc[100];   //[Npvcscnd]
   Int_t           cluster_ncand;
   Float_t         cluster_tstart[10];   //[cluster_ncand]
   Float_t         cluster_tend[10];   //[cluster_ncand]
   Int_t           cluster_nhits[10];   //[cluster_ncand]
   Float_t         cluster_totq[10];   //[cluster_ncand]
   Int_t           cluster_goodflag[10];   //[cluster_ncand]
   Int_t           cluster_npeaks[10][6];   //[cluster_ncand]
   Int_t           cluster_ipeak[10][6][10];   //[cluster_ncand]
   Float_t         cluster_timeofpeak[10][6][10];   //[cluster_ncand]
   Int_t           muechk_ncand[6];
   Float_t         muechk_toya[6];
   Float_t         muechk_tpeak[6][10];
   Int_t           muechk_bg[6][10];
   Float_t         muechk_mean[6][10];
   Float_t         muechk_excess[6][10];
   Float_t         muechk_signif[6][10];
   Int_t           muechk_icluster[6][10];
   Float_t         trgoff;
   Int_t           fqntwnd;
   Int_t           fqtwnd_iclstr[10];   //[fqntwnd]
   Int_t           fqtwnd_npeak[10];   //[fqntwnd]
   Float_t         fqtwnd_prftt0[10];   //[fqntwnd]
   Float_t         fqtwnd_prftpos[10][3];   //[fqntwnd]
   Float_t         fqtwnd[10][2];   //[fqntwnd]
   Float_t         fqtwnd_peakt0[10][10];   //[fqntwnd]
   Float_t         fqtwnd_peakiness[10][10];   //[fqntwnd]
   Int_t           fqnse;
   Int_t           fqitwnd[10];   //[fqnse]
   Int_t           fqipeak[10];   //[fqnse]
   Int_t           fqnhitpmt[10];   //[fqnse]
   Float_t         fqtotq[10];   //[fqnse]
   Float_t         fq0rtotmu[10];   //[fqnse]
   Float_t         fq0rnll[10];   //[fqnse]
   Int_t           fqn50[10];   //[fqnse]
   Float_t         fqq50[10];   //[fqnse]
   Int_t           fq1rpcflg[10][7];   //[fqnse]
   Float_t         fq1rmom[10][7];   //[fqnse]
   Float_t         fq1rt0[10][7];   //[fqnse]
   Float_t         fq1rtotmu[10][7];   //[fqnse]
   Float_t         fq1rnll[10][7];   //[fqnse]
   Float_t         fq1rpos[10][7][3];   //[fqnse]
   Float_t         fq1rdir[10][7][3];   //[fqnse]
   Float_t         fq1rdconv[10][7];   //[fqnse]
   Float_t         fq1reloss[10][7];   //[fqnse]
   Int_t           fqpi0pcflg[2];
   Float_t         fqpi0mom1[2];
   Float_t         fqpi0mom2[2];
   Float_t         fqpi0momtot[2];
   Float_t         fqpi0dconv1[2];
   Float_t         fqpi0dconv2[2];
   Float_t         fqpi0t0[2];
   Float_t         fqpi0totmu[2];
   Float_t         fqpi0nll[2];
   Float_t         fqpi0mass[2];
   Float_t         fqpi0photangle[2];
   Float_t         fqpi0pos[2][3];
   Float_t         fqpi0dir1[2][3];
   Float_t         fqpi0dir2[2][3];
   Float_t         fqpi0dirtot[2][3];
   Int_t           fqpmgpcflg[2];
   Float_t         fqpmgmom1[2];
   Float_t         fqpmgmom2[2];
   Float_t         fqpmgt01[2];
   Float_t         fqpmgt02[2];
   Float_t         fqpmgtotmu[2];
   Float_t         fqpmgnll[2];
   Float_t         fqpmgpos1[2][3];
   Float_t         fqpmgpos2[2][3];
   Float_t         fqpmgdir1[2][3];
   Float_t         fqpmgdir2[2][3];
   Int_t           fqnmrfit;
   Int_t           fqmrifit[200];   //[fqnmrfit]
   Int_t           fqmrnring[200];   //[fqnmrfit]
   Int_t           fqmrpcflg[200];   //[fqnmrfit]
   Float_t         fqmrnll[200];   //[fqnmrfit]
   Float_t         fqmrtotmu[200];   //[fqnmrfit]
   Int_t           fqmrpid[200][6];   //[fqnmrfit]
   Float_t         fqmrmom[200][6];   //[fqnmrfit]
   Float_t         fqmrdconv[200][6];   //[fqnmrfit]
   Float_t         fqmreloss[200][6];   //[fqnmrfit]
   Float_t         fqmrt0[200][6];   //[fqnmrfit]
   Float_t         fqmrpos[200][6][3];   //[fqnmrfit]
   Float_t         fqmrdir[200][6][3];   //[fqnmrfit]
   Int_t           fqmsnfit;
   Int_t           fqmspcflg[5];   //[fqmsnfit]
   Int_t           fqmsnseg[5];   //[fqmsnfit]
   Int_t           fqmspid[5];   //[fqmsnfit]
   Int_t           fqmsifit[5];   //[fqmsnfit]
   Int_t           fqmsimer[5];   //[fqmsnfit]
   Float_t         fqmstotmu[5];   //[fqmsnfit]
   Float_t         fqmsnll[5];   //[fqmsnfit]
   Float_t         fqmsmom[5][20];   //[fqmsnfit]
   Float_t         fqmseloss[5][20];   //[fqmsnfit]
   Float_t         fqmst0[5][20];   //[fqmsnfit]
   Float_t         fqmspos[5][20][3];   //[fqmsnfit]
   Float_t         fqmsdir[5][20][3];   //[fqmsnfit]
   Int_t           fqtestn1r;
   Int_t           fqtest1rstage[30];   //[fqtestn1r]
   Int_t           fqtest1rse[30];   //[fqtestn1r]
   Int_t           fqtest1rpid[30];   //[fqtestn1r]
   Int_t           fqtest1rpcflg[30];   //[fqtestn1r]
   Float_t         fqtest1rmom[30];   //[fqtestn1r]
   Float_t         fqtest1rt0[30];   //[fqtestn1r]
   Float_t         fqtest1rtotmu[30];   //[fqtestn1r]
   Float_t         fqtest1rnll[30];   //[fqtestn1r]
   Float_t         fqtest1rpos[30][3];   //[fqtestn1r]
   Float_t         fqtest1rdir[30][3];   //[fqtestn1r]
   Float_t         fqtest1rdconv[30];   //[fqtestn1r]
   Float_t         fqtest1reloss[30];   //[fqtestn1r]
   Int_t           fqtestnpi0;
   Int_t           fqtestpi0stage[30];   //[fqtestnpi0]
   Int_t           fqtestpi0pcflg[30];   //[fqtestnpi0]
   Float_t         fqtestpi0mom1[30];   //[fqtestnpi0]
   Float_t         fqtestpi0mom2[30];   //[fqtestnpi0]
   Float_t         fqtestpi0momtot[30];   //[fqtestnpi0]
   Float_t         fqtestpi0dconv1[30];   //[fqtestnpi0]
   Float_t         fqtestpi0dconv2[30];   //[fqtestnpi0]
   Float_t         fqtestpi0t0[30];   //[fqtestnpi0]
   Float_t         fqtestpi0totmu[30];   //[fqtestnpi0]
   Float_t         fqtestpi0nll[30];   //[fqtestnpi0]
   Float_t         fqtestpi0mass[30];   //[fqtestnpi0]
   Float_t         fqtestpi0photangle[30];   //[fqtestnpi0]
   Float_t         fqtestpi0pos[30][3];   //[fqtestnpi0]
   Float_t         fqtestpi0dir1[30][3];   //[fqtestnpi0]
   Float_t         fqtestpi0dir2[30][3];   //[fqtestnpi0]
   Float_t         fqtestpi0dirtot[30][3];   //[fqtestnpi0]
   Float_t         prid_ang;
   Float_t         prid_len;
   Float_t         prid_pmom;
   Float_t         prid_mang;
   Float_t         prid_plike;
   Float_t         prid_mlike;
   Int_t           prid_nused;
   Float_t         dirtotepi[3];
   Float_t         dirtotenpi[3];
   Float_t         dirtotmue[3];
   Float_t         dirsum[3];
   Float_t         etot;
   Float_t         etotepi;
   Float_t         etotenpi;
   Float_t         etotmue;
   Float_t         oscweight[2][4];
   Float_t         oscwgt;
   Float_t         live;
   Float_t         sacth;
   Float_t         sactg;
   Float_t         sacth1d;
   Float_t         ent_pos[3];
   Float_t         ent_dir[3];
   Float_t         length;
   Float_t         tr_mom1;
   Float_t         A_ent_mom;
   Float_t         A_ent_pos[3];
   Float_t         A_ent_dir[3];
   Float_t         A_ext_mom;
   Float_t         A_ext_pos[3];
   Float_t         A_ext_dir[3];
   Float_t         Fit_pos[3];
   Float_t         Fit_dir[3];
   Float_t         Fit_len;
   Float_t         Fit_mom;
   Int_t           Fit_pid;
   Int_t           Um_ehit8m;
   Int_t           Um_ohit8m;
   Float_t         Um_qent;
   Float_t         Sh_chi1p;
   Float_t         Sh_delta;
   Float_t         Sh_mean;
   Float_t         Sh_meanq;
   Int_t           Sh_id;
   Float_t         Tr_stop[3];
   Float_t         Tr_mom;
   Float_t         Tr_len;
   Float_t         Tr_len1;
   Int_t           Pid_flg;
   Float_t         Crs1;
   Float_t         Crs2;
   Int_t           iclass;
   Int_t           mu_class;
   Int_t           mu_dec;
   Float_t         mu_dir[3];
   Float_t         mu_pos[3];
   Float_t         mu_good;
   Int_t           history;
   Int_t           Pdst;
   Int_t           idoff;
   Float_t         anthit;
   Int_t           idseq;
   Float_t         tstfrac;
   Int_t           judge;
   Float_t         Upcrs1;
   Float_t         Upcrs2;
   Float_t         lst;
   Int_t           jd;
   Float_t         fjd;
   Float_t         alt;
   Float_t         azi;
   Float_t         ra;
   Float_t         dec;
   Float_t         glat;
   Float_t         glong;

   // List of branches
   TBranch        *b_nring;   //!
   TBranch        *b_nrun;   //!
   TBranch        *b_nev;   //!
   TBranch        *b_nsub;   //!
   TBranch        *b_cate;   //!
   TBranch        *b_potot;   //!
   TBranch        *b_nhit;   //!
   TBranch        *b_pomax;   //!
   TBranch        *b_potota;   //!
   TBranch        *b_nhita;   //!
   TBranch        *b_nhitac;   //!
   TBranch        *b_pomaxa;   //!
   TBranch        *b_wall;   //!
   TBranch        *b_evis;   //!
   TBranch        *b_rtsum;   //!
   TBranch        *b_rtmax;   //!
   TBranch        *b_wlen;   //!
   TBranch        *b_ip;   //!
   TBranch        *b_pos;   //!
   TBranch        *b_dir;   //!
   TBranch        *b_dirtot;   //!
   TBranch        *b_ang;   //!
   TBranch        *b_rtot;   //!
   TBranch        *b_amom;   //!
   TBranch        *b_rtote;   //!
   TBranch        *b_amome;   //!
   TBranch        *b_rtotm;   //!
   TBranch        *b_amomm;   //!
   TBranch        *b_nsube;   //!
   TBranch        *b_ndcy;   //!
   TBranch        *b_ngate;   //!
   TBranch        *b_nbye;   //!
   TBranch        *b_probms;   //!
   TBranch        *b_prmslg;   //!
   TBranch        *b_date;   //!
   TBranch        *b_time;   //!
   TBranch        *b_elpsday;   //!
   TBranch        *b_numpo;   //!
   TBranch        *b_apos;   //!
   TBranch        *b_adir;   //!
   TBranch        *b_aang;   //!
   TBranch        *b_agood;   //!
   TBranch        *b_wgain;   //!
   TBranch        *b_nbad;   //!
   TBranch        *b_nbada;   //!
   TBranch        *b_msdir;   //!
   TBranch        *b_amomp;   //!
   TBranch        *b_ange;   //!
   TBranch        *b_angm;   //!
   TBranch        *b_angp;   //!
   TBranch        *b_ntot;   //!
   TBranch        *b_probth;   //!
   TBranch        *b_probpt;   //!
   TBranch        *b_pi0like;   //!
   TBranch        *b_pi0_e;   //!
   TBranch        *b_pi0_dir;   //!
   TBranch        *b_pi0mass;   //!
   TBranch        *b_evisold;   //!
   TBranch        *b_evisoldxe;   //!
   TBranch        *b_evisnew;   //!
   TBranch        *b_nmue;   //!
   TBranch        *b_etype;   //!
   TBranch        *b_etime;   //!
   TBranch        *b_epos;   //!
   TBranch        *b_edir;   //!
   TBranch        *b_egood;   //!
   TBranch        *b_ehit;   //!
   TBranch        *b_mueprob;   //!
   TBranch        *b_nmue_sel;   //!
   TBranch        *b_Rnring;   //!
   TBranch        *b_Rdir;   //!
   TBranch        *b_Rang;   //!
   TBranch        *b_Riring;   //!
   TBranch        *b_Rtwout;   //!
   TBranch        *b_Rtwith;   //!
   TBranch        *b_Alwout;   //!
   TBranch        *b_Alwith;   //!
   TBranch        *b_Qsmi;   //!
   TBranch        *b_Qsmo;   //!
   TBranch        *b_Qexi;   //!
   TBranch        *b_Qexo;   //!
   TBranch        *b_Pe5d;   //!
   TBranch        *b_En5d;   //!
   TBranch        *b_Eh5d;   //!
   TBranch        *b_Pe5do;   //!
   TBranch        *b_En5do;   //!
   TBranch        *b_Eh5do;   //!
   TBranch        *b_Rtadd;   //!
   TBranch        *b_Pdgeta;   //!
   TBranch        *b_Pd5d;   //!
   TBranch        *b_Pdthre;   //!
   TBranch        *b_Pd5do;   //!
   TBranch        *b_Delpd;   //!
   TBranch        *b_Ropena;   //!
   TBranch        *b_Maxth;   //!
   TBranch        *b_Pkang;   //!
   TBranch        *b_Qrfct;   //!
   TBranch        *b_Pdfct;   //!
   TBranch        *b_Pkfct;   //!
   TBranch        *b_Agfct;   //!
   TBranch        *b_Dlfct;   //!
   TBranch        *b_Iflag;   //!
   TBranch        *b_Pmfct;   //!
   TBranch        *b_Imfct;   //!
   TBranch        *b_Rilike;   //!
   TBranch        *b_ri_ver;   //!
   TBranch        *b_ri_pid;   //!
   TBranch        *b_ri_nring;   //!
   TBranch        *b_ri_flag;   //!
   TBranch        *b_ri_dlfct;   //!
   TBranch        *b_ri_pdfct;   //!
   TBranch        *b_ri_pkfct;   //!
   TBranch        *b_ri_vafct;   //!
   TBranch        *b_ri_total;   //!
   TBranch        *b_ri_dir;   //!
   TBranch        *b_ri_imfct;   //!
   TBranch        *b_ri_pmfct;   //!
   TBranch        *b_npar;   //!
   TBranch        *b_wallv;   //!
   TBranch        *b_ipv;   //!
   TBranch        *b_posv;   //!
   TBranch        *b_dirv;   //!
   TBranch        *b_pmomv;   //!
   TBranch        *b_light_flag;   //!
   TBranch        *b_npar2;   //!
   TBranch        *b_wallv2;   //!
   TBranch        *b_ipv2;   //!
   TBranch        *b_iorg;   //!
   TBranch        *b_posv2;   //!
   TBranch        *b_dirv2;   //!
   TBranch        *b_pmomv2;   //!
   TBranch        *b_numnu;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_ipnu;   //!
   TBranch        *b_pnu;   //!
   TBranch        *b_dirnu;   //!
   TBranch        *b_flxg;   //!
   TBranch        *b_flxh01;   //!
   TBranch        *b_kflux;   //!
   TBranch        *b_bs71;   //!
   TBranch        *b_bs74;   //!
   TBranch        *b_flxf;   //!
   TBranch        *b_flxh1d;   //!
   TBranch        *b_flxb03;   //!
   TBranch        *b_flxf03;   //!
   TBranch        *b_flxh06;   //!
   TBranch        *b_flxh11;   //!
   TBranch        *b_scan;   //!
   TBranch        *b_nuceff_version;   //!
   TBranch        *b_charge_exchange;   //!
   TBranch        *b_absorbed;   //!
   TBranch        *b_multipi_gen;   //!
   TBranch        *b_scattering;   //!
   TBranch        *b_piless_dcy;   //!
   TBranch        *b_nscndprt;   //!
   TBranch        *b_itrkscnd;   //!
   TBranch        *b_istakscnd;   //!
   TBranch        *b_vtxscnd;   //!
   TBranch        *b_pscnd;   //!
   TBranch        *b_iprtscnd;   //!
   TBranch        *b_tscnd;   //!
   TBranch        *b_iprntprt;   //!
   TBranch        *b_lmecscnd;   //!
   TBranch        *b_iprnttrk;   //!
   TBranch        *b_iorgprt;   //!
   TBranch        *b_pprnt;   //!
   TBranch        *b_iflgscnd;   //!
   TBranch        *b_pprntinit;   //!
   TBranch        *b_vtxprnt;   //!
   TBranch        *b_iprntidx;   //!
   TBranch        *b_nchilds;   //!
   TBranch        *b_ichildidx;   //!
   TBranch        *b_Npvcscnd;   //!
   TBranch        *b_Nchildsvc;   //!
   TBranch        *b_Ichildidxvc;   //!
   TBranch        *b_cluster_ncand;   //!
   TBranch        *b_cluster_tstart;   //!
   TBranch        *b_cluster_tend;   //!
   TBranch        *b_cluster_nhits;   //!
   TBranch        *b_cluster_totq;   //!
   TBranch        *b_cluster_goodflag;   //!
   TBranch        *b_cluster_npeaks;   //!
   TBranch        *b_cluster_ipeak;   //!
   TBranch        *b_cluster_timeofpeak;   //!
   TBranch        *b_muechk_ncand;   //!
   TBranch        *b_muechk_toya;   //!
   TBranch        *b_muechk_tpeak;   //!
   TBranch        *b_muechk_bg;   //!
   TBranch        *b_muechk_mean;   //!
   TBranch        *b_muechk_excess;   //!
   TBranch        *b_muechk_signif;   //!
   TBranch        *b_muechk_icluster;   //!
   TBranch        *b_trgoff;   //!
   TBranch        *b_fqntwnd;   //!
   TBranch        *b_fqtwnd_iclstr;   //!
   TBranch        *b_fqtwnd_npeak;   //!
   TBranch        *b_fqtwnd_prftt0;   //!
   TBranch        *b_fqtwnd_prftpos;   //!
   TBranch        *b_fqtwnd;   //!
   TBranch        *b_fqtwnd_peakt0;   //!
   TBranch        *b_fqtwnd_peakiness;   //!
   TBranch        *b_fqnse;   //!
   TBranch        *b_fqitwnd;   //!
   TBranch        *b_fqipeak;   //!
   TBranch        *b_fqnhitpmt;   //!
   TBranch        *b_fqtotq;   //!
   TBranch        *b_fq0rtotmu;   //!
   TBranch        *b_fq0rnll;   //!
   TBranch        *b_fqn50;   //!
   TBranch        *b_fqq50;   //!
   TBranch        *b_fq1rpcflg;   //!
   TBranch        *b_fq1rmom;   //!
   TBranch        *b_fq1rt0;   //!
   TBranch        *b_fq1rtotmu;   //!
   TBranch        *b_fq1rnll;   //!
   TBranch        *b_fq1rpos;   //!
   TBranch        *b_fq1rdir;   //!
   TBranch        *b_fq1rdconv;   //!
   TBranch        *b_fq1reloss;   //!
   TBranch        *b_fqpi0pcflg;   //!
   TBranch        *b_fqpi0mom1;   //!
   TBranch        *b_fqpi0mom2;   //!
   TBranch        *b_fqpi0momtot;   //!
   TBranch        *b_fqpi0dconv1;   //!
   TBranch        *b_fqpi0dconv2;   //!
   TBranch        *b_fqpi0t0;   //!
   TBranch        *b_fqpi0totmu;   //!
   TBranch        *b_fqpi0nll;   //!
   TBranch        *b_fqpi0mass;   //!
   TBranch        *b_fqpi0photangle;   //!
   TBranch        *b_fqpi0pos;   //!
   TBranch        *b_fqpi0dir1;   //!
   TBranch        *b_fqpi0dir2;   //!
   TBranch        *b_fqpi0dirtot;   //!
   TBranch        *b_fqpmgpcflg;   //!
   TBranch        *b_fqpmgmom1;   //!
   TBranch        *b_fqpmgmom2;   //!
   TBranch        *b_fqpmgt01;   //!
   TBranch        *b_fqpmgt02;   //!
   TBranch        *b_fqpmgtotmu;   //!
   TBranch        *b_fqpmgnll;   //!
   TBranch        *b_fqpmgpos1;   //!
   TBranch        *b_fqpmgpos2;   //!
   TBranch        *b_fqpmgdir1;   //!
   TBranch        *b_fqpmgdir2;   //!
   TBranch        *b_fqnmrfit;   //!
   TBranch        *b_fqmrifit;   //!
   TBranch        *b_fqmrnring;   //!
   TBranch        *b_fqmrpcflg;   //!
   TBranch        *b_fqmrnll;   //!
   TBranch        *b_fqmrtotmu;   //!
   TBranch        *b_fqmrpid;   //!
   TBranch        *b_fqmrmom;   //!
   TBranch        *b_fqmrdconv;   //!
   TBranch        *b_fqmreloss;   //!
   TBranch        *b_fqmrt0;   //!
   TBranch        *b_fqmrpos;   //!
   TBranch        *b_fqmrdir;   //!
   TBranch        *b_fqmsnfit;   //!
   TBranch        *b_fqmspcflg;   //!
   TBranch        *b_fqmsnseg;   //!
   TBranch        *b_fqmspid;   //!
   TBranch        *b_fqmsifit;   //!
   TBranch        *b_fqmsimer;   //!
   TBranch        *b_fqmstotmu;   //!
   TBranch        *b_fqmsnll;   //!
   TBranch        *b_fqmsmom;   //!
   TBranch        *b_fqmseloss;   //!
   TBranch        *b_fqmst0;   //!
   TBranch        *b_fqmspos;   //!
   TBranch        *b_fqmsdir;   //!
   TBranch        *b_fqtestn1r;   //!
   TBranch        *b_fqtest1rstage;   //!
   TBranch        *b_fqtest1rse;   //!
   TBranch        *b_fqtest1rpid;   //!
   TBranch        *b_fqtest1rpcflg;   //!
   TBranch        *b_fqtest1rmom;   //!
   TBranch        *b_fqtest1rt0;   //!
   TBranch        *b_fqtest1rtotmu;   //!
   TBranch        *b_fqtest1rnll;   //!
   TBranch        *b_fqtest1rpos;   //!
   TBranch        *b_fqtest1rdir;   //!
   TBranch        *b_fqtest1rdconv;   //!
   TBranch        *b_fqtest1reloss;   //!
   TBranch        *b_fqtestnpi0;   //!
   TBranch        *b_fqtestpi0stage;   //!
   TBranch        *b_fqtestpi0pcflg;   //!
   TBranch        *b_fqtestpi0mom1;   //!
   TBranch        *b_fqtestpi0mom2;   //!
   TBranch        *b_fqtestpi0momtot;   //!
   TBranch        *b_fqtestpi0dconv1;   //!
   TBranch        *b_fqtestpi0dconv2;   //!
   TBranch        *b_fqtestpi0t0;   //!
   TBranch        *b_fqtestpi0totmu;   //!
   TBranch        *b_fqtestpi0nll;   //!
   TBranch        *b_fqtestpi0mass;   //!
   TBranch        *b_fqtestpi0photangle;   //!
   TBranch        *b_fqtestpi0pos;   //!
   TBranch        *b_fqtestpi0dir1;   //!
   TBranch        *b_fqtestpi0dir2;   //!
   TBranch        *b_fqtestpi0dirtot;   //!
   TBranch        *b_prid_ang;   //!
   TBranch        *b_prid_len;   //!
   TBranch        *b_prid_pmom;   //!
   TBranch        *b_prid_mang;   //!
   TBranch        *b_prid_plike;   //!
   TBranch        *b_prid_mlike;   //!
   TBranch        *b_prid_nused;   //!
   TBranch        *b_dirtotepi;   //!
   TBranch        *b_dirtotenpi;   //!
   TBranch        *b_dirtotmue;   //!
   TBranch        *b_dirsum;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_etotepi;   //!
   TBranch        *b_etotenpi;   //!
   TBranch        *b_etotmue;   //!
   TBranch        *b_oscweight;   //!
   TBranch        *b_oscwgt;   //!
   TBranch        *b_live;   //!
   TBranch        *b_sacth;   //!
   TBranch        *b_sactg;   //!
   TBranch        *b_sacth1d;   //!
   TBranch        *b_ent_pos;   //!
   TBranch        *b_ent_dir;   //!
   TBranch        *b_length;   //!
   TBranch        *b_tr_mom1;   //!
   TBranch        *b_A_ent_mom;   //!
   TBranch        *b_A_ent_pos;   //!
   TBranch        *b_A_ent_dir;   //!
   TBranch        *b_A_ext_mom;   //!
   TBranch        *b_A_ext_pos;   //!
   TBranch        *b_A_ext_dir;   //!
   TBranch        *b_Fit_pos;   //!
   TBranch        *b_Fit_dir;   //!
   TBranch        *b_Fit_len;   //!
   TBranch        *b_Fit_mom;   //!
   TBranch        *b_Fit_pid;   //!
   TBranch        *b_Um_ehit8m;   //!
   TBranch        *b_Um_ohit8m;   //!
   TBranch        *b_Um_qent;   //!
   TBranch        *b_Sh_chi1p;   //!
   TBranch        *b_Sh_delta;   //!
   TBranch        *b_Sh_mean;   //!
   TBranch        *b_Sh_meanq;   //!
   TBranch        *b_Sh_id;   //!
   TBranch        *b_Tr_stop;   //!
   TBranch        *b_Tr_mom;   //!
   TBranch        *b_Tr_len;   //!
   TBranch        *b_Tr_len1;   //!
   TBranch        *b_Pid_flg;   //!
   TBranch        *b_Crs1;   //!
   TBranch        *b_Crs2;   //!
   TBranch        *b_iclass;   //!
   TBranch        *b_mu_class;   //!
   TBranch        *b_mu_dec;   //!
   TBranch        *b_mu_dir;   //!
   TBranch        *b_mu_pos;   //!
   TBranch        *b_mu_good;   //!
   TBranch        *b_history;   //!
   TBranch        *b_Pdst;   //!
   TBranch        *b_idoff;   //!
   TBranch        *b_anthit;   //!
   TBranch        *b_idseq;   //!
   TBranch        *b_tstfrac;   //!
   TBranch        *b_judge;   //!
   TBranch        *b_Upcrs1;   //!
   TBranch        *b_Upcrs2;   //!
   TBranch        *b_lst;   //!
   TBranch        *b_jd;   //!
   TBranch        *b_fjd;   //!
   TBranch        *b_alt;   //!
   TBranch        *b_azi;   //!
   TBranch        *b_ra;   //!
   TBranch        *b_dec;   //!
   TBranch        *b_glat;   //!
   TBranch        *b_glong;   //!

   NewVars(TString inputfile="/Users/santucci/Dropbox/PhD/SK/fiTQun_analysis/Knu_muGamma/atmnu/multiring/files/pdk/pdk_100k.root", TTree *tree=0);
   virtual ~NewVars();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString outputfile, int option);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NewVars_cxx
NewVars::NewVars(TString inputfile, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("h1",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("h1","");
      chain->Add(inputfile);
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

NewVars::~NewVars()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NewVars::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NewVars::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NewVars::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nring", &nring, &b_nring);
   fChain->SetBranchAddress("nrun", &nrun, &b_nrun);
   fChain->SetBranchAddress("nev", &nev, &b_nev);
   fChain->SetBranchAddress("nsub", &nsub, &b_nsub);
   fChain->SetBranchAddress("cate", &cate, &b_cate);
   fChain->SetBranchAddress("potot", &potot, &b_potot);
   fChain->SetBranchAddress("nhit", &nhit, &b_nhit);
   fChain->SetBranchAddress("pomax", &pomax, &b_pomax);
   fChain->SetBranchAddress("potota", &potota, &b_potota);
   fChain->SetBranchAddress("nhita", &nhita, &b_nhita);
   fChain->SetBranchAddress("nhitac", &nhitac, &b_nhitac);
   fChain->SetBranchAddress("pomaxa", &pomaxa, &b_pomaxa);
   fChain->SetBranchAddress("wall", &wall, &b_wall);
   fChain->SetBranchAddress("evis", &evis, &b_evis);
   fChain->SetBranchAddress("rtsum", &rtsum, &b_rtsum);
   fChain->SetBranchAddress("rtmax", &rtmax, &b_rtmax);
   fChain->SetBranchAddress("wlen", &wlen, &b_wlen);
   fChain->SetBranchAddress("ip", ip, &b_ip);
   fChain->SetBranchAddress("pos", pos, &b_pos);
   fChain->SetBranchAddress("dir", dir, &b_dir);
   fChain->SetBranchAddress("dirtot", dirtot, &b_dirtot);
   fChain->SetBranchAddress("ang", ang, &b_ang);
   fChain->SetBranchAddress("rtot", rtot, &b_rtot);
   fChain->SetBranchAddress("amom", amom, &b_amom);
   fChain->SetBranchAddress("rtote", rtote, &b_rtote);
   fChain->SetBranchAddress("amome", amome, &b_amome);
   fChain->SetBranchAddress("rtotm", rtotm, &b_rtotm);
   fChain->SetBranchAddress("amomm", amomm, &b_amomm);
   fChain->SetBranchAddress("nsube", &nsube, &b_nsube);
   fChain->SetBranchAddress("ndcy", &ndcy, &b_ndcy);
   fChain->SetBranchAddress("ngate", &ngate, &b_ngate);
   fChain->SetBranchAddress("nbye", &nbye, &b_nbye);
   fChain->SetBranchAddress("probms", probms, &b_probms);
   fChain->SetBranchAddress("prmslg", prmslg, &b_prmslg);
   fChain->SetBranchAddress("date", date, &b_date);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("elpsday", &elpsday, &b_elpsday);
   fChain->SetBranchAddress("numpo", numpo, &b_numpo);
   fChain->SetBranchAddress("apos", apos, &b_apos);
   fChain->SetBranchAddress("adir", adir, &b_adir);
   fChain->SetBranchAddress("aang", &aang, &b_aang);
   fChain->SetBranchAddress("agood", &agood, &b_agood);
   fChain->SetBranchAddress("wgain", &wgain, &b_wgain);
   fChain->SetBranchAddress("nbad", &nbad, &b_nbad);
   fChain->SetBranchAddress("nbada", &nbada, &b_nbada);
   fChain->SetBranchAddress("msdir", msdir, &b_msdir);
   fChain->SetBranchAddress("amomp", amomp, &b_amomp);
   fChain->SetBranchAddress("ange", ange, &b_ange);
   fChain->SetBranchAddress("angm", angm, &b_angm);
   fChain->SetBranchAddress("angp", angp, &b_angp);
   fChain->SetBranchAddress("ntot", ntot, &b_ntot);
   fChain->SetBranchAddress("probth", probth, &b_probth);
   fChain->SetBranchAddress("probpt", probpt, &b_probpt);
   fChain->SetBranchAddress("pi0like", pi0like, &b_pi0like);
   fChain->SetBranchAddress("pi0_e", pi0_e, &b_pi0_e);
   fChain->SetBranchAddress("pi0_dir", pi0_dir, &b_pi0_dir);
   fChain->SetBranchAddress("pi0mass", pi0mass, &b_pi0mass);
   fChain->SetBranchAddress("evisold", &evisold, &b_evisold);
   fChain->SetBranchAddress("evisoldxe", &evisoldxe, &b_evisoldxe);
   fChain->SetBranchAddress("evisnew", &evisnew, &b_evisnew);
   fChain->SetBranchAddress("nmue", &nmue, &b_nmue);
   fChain->SetBranchAddress("etype", etype, &b_etype);
   fChain->SetBranchAddress("etime", etime, &b_etime);
   fChain->SetBranchAddress("epos", epos, &b_epos);
   fChain->SetBranchAddress("edir", edir, &b_edir);
   fChain->SetBranchAddress("egood", egood, &b_egood);
   fChain->SetBranchAddress("ehit", ehit, &b_ehit);
   fChain->SetBranchAddress("mueprob", mueprob, &b_mueprob);
   fChain->SetBranchAddress("nmue_sel", &nmue_sel, &b_nmue_sel);
   fChain->SetBranchAddress("Rnring", &Rnring, &b_Rnring);
   fChain->SetBranchAddress("Rdir", Rdir, &b_Rdir);
   fChain->SetBranchAddress("Rang", Rang, &b_Rang);
   fChain->SetBranchAddress("Riring", &Riring, &b_Riring);
   fChain->SetBranchAddress("Rtwout", Rtwout, &b_Rtwout);
   fChain->SetBranchAddress("Rtwith", Rtwith, &b_Rtwith);
   fChain->SetBranchAddress("Alwout", &Alwout, &b_Alwout);
   fChain->SetBranchAddress("Alwith", &Alwith, &b_Alwith);
   fChain->SetBranchAddress("Qsmi", &Qsmi, &b_Qsmi);
   fChain->SetBranchAddress("Qsmo", &Qsmo, &b_Qsmo);
   fChain->SetBranchAddress("Qexi", &Qexi, &b_Qexi);
   fChain->SetBranchAddress("Qexo", &Qexo, &b_Qexo);
   fChain->SetBranchAddress("Pe5d", &Pe5d, &b_Pe5d);
   fChain->SetBranchAddress("En5d", &En5d, &b_En5d);
   fChain->SetBranchAddress("Eh5d", &Eh5d, &b_Eh5d);
   fChain->SetBranchAddress("Pe5do", &Pe5do, &b_Pe5do);
   fChain->SetBranchAddress("En5do", &En5do, &b_En5do);
   fChain->SetBranchAddress("Eh5do", &Eh5do, &b_Eh5do);
   fChain->SetBranchAddress("Rtadd", &Rtadd, &b_Rtadd);
   fChain->SetBranchAddress("Pdgeta", &Pdgeta, &b_Pdgeta);
   fChain->SetBranchAddress("Pd5d", &Pd5d, &b_Pd5d);
   fChain->SetBranchAddress("Pdthre", &Pdthre, &b_Pdthre);
   fChain->SetBranchAddress("Pd5do", &Pd5do, &b_Pd5do);
   fChain->SetBranchAddress("Delpd", &Delpd, &b_Delpd);
   fChain->SetBranchAddress("Ropena", Ropena, &b_Ropena);
   fChain->SetBranchAddress("Maxth", &Maxth, &b_Maxth);
   fChain->SetBranchAddress("Pkang", &Pkang, &b_Pkang);
   fChain->SetBranchAddress("Qrfct", &Qrfct, &b_Qrfct);
   fChain->SetBranchAddress("Pdfct", &Pdfct, &b_Pdfct);
   fChain->SetBranchAddress("Pkfct", &Pkfct, &b_Pkfct);
   fChain->SetBranchAddress("Agfct", &Agfct, &b_Agfct);
   fChain->SetBranchAddress("Dlfct", &Dlfct, &b_Dlfct);
   fChain->SetBranchAddress("Iflag", &Iflag, &b_Iflag);
   fChain->SetBranchAddress("Pmfct", &Pmfct, &b_Pmfct);
   fChain->SetBranchAddress("Imfct", &Imfct, &b_Imfct);
   fChain->SetBranchAddress("Rilike", &Rilike, &b_Rilike);
   fChain->SetBranchAddress("ri_ver", &ri_ver, &b_ri_ver);
   fChain->SetBranchAddress("ri_pid", &ri_pid, &b_ri_pid);
   fChain->SetBranchAddress("ri_nring", &ri_nring, &b_ri_nring);
   fChain->SetBranchAddress("ri_flag", ri_flag, &b_ri_flag);
   fChain->SetBranchAddress("ri_dlfct", ri_dlfct, &b_ri_dlfct);
   fChain->SetBranchAddress("ri_pdfct", ri_pdfct, &b_ri_pdfct);
   fChain->SetBranchAddress("ri_pkfct", ri_pkfct, &b_ri_pkfct);
   fChain->SetBranchAddress("ri_vafct", ri_vafct, &b_ri_vafct);
   fChain->SetBranchAddress("ri_total", ri_total, &b_ri_total);
   fChain->SetBranchAddress("ri_dir", ri_dir, &b_ri_dir);
   fChain->SetBranchAddress("ri_imfct", ri_imfct, &b_ri_imfct);
   fChain->SetBranchAddress("ri_pmfct", ri_pmfct, &b_ri_pmfct);
   fChain->SetBranchAddress("npar", &npar, &b_npar);
   fChain->SetBranchAddress("wallv", &wallv, &b_wallv);
   fChain->SetBranchAddress("ipv", ipv, &b_ipv);
   fChain->SetBranchAddress("posv", posv, &b_posv);
   fChain->SetBranchAddress("dirv", dirv, &b_dirv);
   fChain->SetBranchAddress("pmomv", pmomv, &b_pmomv);
   fChain->SetBranchAddress("light_flag", light_flag, &b_light_flag);
   fChain->SetBranchAddress("npar2", &npar2, &b_npar2);
   fChain->SetBranchAddress("wallv2", wallv2, &b_wallv2);
   fChain->SetBranchAddress("ipv2", ipv2, &b_ipv2);
   fChain->SetBranchAddress("iorg", iorg, &b_iorg);
   fChain->SetBranchAddress("posv2", posv2, &b_posv2);
   fChain->SetBranchAddress("dirv2", dirv2, &b_dirv2);
   fChain->SetBranchAddress("pmomv2", pmomv2, &b_pmomv2);
   fChain->SetBranchAddress("numnu", &numnu, &b_numnu);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("ipnu", ipnu, &b_ipnu);
   fChain->SetBranchAddress("pnu", pnu, &b_pnu);
   fChain->SetBranchAddress("dirnu", dirnu, &b_dirnu);
   fChain->SetBranchAddress("flxg", flxg, &b_flxg);
   fChain->SetBranchAddress("flxh01", flxh01, &b_flxh01);
   fChain->SetBranchAddress("kflux", kflux, &b_kflux);
   fChain->SetBranchAddress("bs71", bs71, &b_bs71);
   fChain->SetBranchAddress("bs74", bs74, &b_bs74);
   fChain->SetBranchAddress("flxf", flxf, &b_flxf);
   fChain->SetBranchAddress("flxh1d", flxh1d, &b_flxh1d);
   fChain->SetBranchAddress("flxb03", flxb03, &b_flxb03);
   fChain->SetBranchAddress("flxf03", flxf03, &b_flxf03);
   fChain->SetBranchAddress("flxh06", flxh06, &b_flxh06);
   fChain->SetBranchAddress("flxh11", flxh11, &b_flxh11);
   fChain->SetBranchAddress("scan", scan, &b_scan);
   fChain->SetBranchAddress("nuceff_version", &nuceff_version, &b_nuceff_version);
   fChain->SetBranchAddress("charge_exchange", &charge_exchange, &b_charge_exchange);
   fChain->SetBranchAddress("absorbed", &absorbed, &b_absorbed);
   fChain->SetBranchAddress("multipi_gen", &multipi_gen, &b_multipi_gen);
   fChain->SetBranchAddress("scattering", &scattering, &b_scattering);
   fChain->SetBranchAddress("piless_dcy", &piless_dcy, &b_piless_dcy);
   fChain->SetBranchAddress("nscndprt", &nscndprt, &b_nscndprt);
   fChain->SetBranchAddress("itrkscnd", itrkscnd, &b_itrkscnd);
   fChain->SetBranchAddress("istakscnd", istakscnd, &b_istakscnd);
   fChain->SetBranchAddress("vtxscnd", vtxscnd, &b_vtxscnd);
   fChain->SetBranchAddress("pscnd", pscnd, &b_pscnd);
   fChain->SetBranchAddress("iprtscnd", iprtscnd, &b_iprtscnd);
   fChain->SetBranchAddress("tscnd", tscnd, &b_tscnd);
   fChain->SetBranchAddress("iprntprt", iprntprt, &b_iprntprt);
   fChain->SetBranchAddress("lmecscnd", lmecscnd, &b_lmecscnd);
   fChain->SetBranchAddress("iprnttrk", iprnttrk, &b_iprnttrk);
   fChain->SetBranchAddress("iorgprt", iorgprt, &b_iorgprt);
   fChain->SetBranchAddress("pprnt", pprnt, &b_pprnt);
   fChain->SetBranchAddress("iflgscnd", iflgscnd, &b_iflgscnd);
   fChain->SetBranchAddress("pprntinit", pprntinit, &b_pprntinit);
   fChain->SetBranchAddress("vtxprnt", vtxprnt, &b_vtxprnt);
   fChain->SetBranchAddress("iprntidx", iprntidx, &b_iprntidx);
   fChain->SetBranchAddress("nchilds", nchilds, &b_nchilds);
   fChain->SetBranchAddress("ichildidx", ichildidx, &b_ichildidx);
   fChain->SetBranchAddress("Npvcscnd", &Npvcscnd, &b_Npvcscnd);
   fChain->SetBranchAddress("Nchildsvc", Nchildsvc, &b_Nchildsvc);
   fChain->SetBranchAddress("Ichildidxvc", Ichildidxvc, &b_Ichildidxvc);
   fChain->SetBranchAddress("cluster_ncand", &cluster_ncand, &b_cluster_ncand);
   fChain->SetBranchAddress("cluster_tstart", cluster_tstart, &b_cluster_tstart);
   fChain->SetBranchAddress("cluster_tend", cluster_tend, &b_cluster_tend);
   fChain->SetBranchAddress("cluster_nhits", cluster_nhits, &b_cluster_nhits);
   fChain->SetBranchAddress("cluster_totq", cluster_totq, &b_cluster_totq);
   fChain->SetBranchAddress("cluster_goodflag", cluster_goodflag, &b_cluster_goodflag);
   fChain->SetBranchAddress("cluster_npeaks", cluster_npeaks, &b_cluster_npeaks);
   fChain->SetBranchAddress("cluster_ipeak", cluster_ipeak, &b_cluster_ipeak);
   fChain->SetBranchAddress("cluster_timeofpeak", cluster_timeofpeak, &b_cluster_timeofpeak);
   fChain->SetBranchAddress("muechk_ncand", muechk_ncand, &b_muechk_ncand);
   fChain->SetBranchAddress("muechk_toya", muechk_toya, &b_muechk_toya);
   fChain->SetBranchAddress("muechk_tpeak", muechk_tpeak, &b_muechk_tpeak);
   fChain->SetBranchAddress("muechk_bg", muechk_bg, &b_muechk_bg);
   fChain->SetBranchAddress("muechk_mean", muechk_mean, &b_muechk_mean);
   fChain->SetBranchAddress("muechk_excess", muechk_excess, &b_muechk_excess);
   fChain->SetBranchAddress("muechk_signif", muechk_signif, &b_muechk_signif);
   fChain->SetBranchAddress("muechk_icluster", muechk_icluster, &b_muechk_icluster);
   fChain->SetBranchAddress("trgoff", &trgoff, &b_trgoff);
   fChain->SetBranchAddress("fqntwnd", &fqntwnd, &b_fqntwnd);
   fChain->SetBranchAddress("fqtwnd_iclstr", fqtwnd_iclstr, &b_fqtwnd_iclstr);
   fChain->SetBranchAddress("fqtwnd_npeak", fqtwnd_npeak, &b_fqtwnd_npeak);
   fChain->SetBranchAddress("fqtwnd_prftt0", fqtwnd_prftt0, &b_fqtwnd_prftt0);
   fChain->SetBranchAddress("fqtwnd_prftpos", fqtwnd_prftpos, &b_fqtwnd_prftpos);
   fChain->SetBranchAddress("fqtwnd", fqtwnd, &b_fqtwnd);
   fChain->SetBranchAddress("fqtwnd_peakt0", fqtwnd_peakt0, &b_fqtwnd_peakt0);
   fChain->SetBranchAddress("fqtwnd_peakiness", fqtwnd_peakiness, &b_fqtwnd_peakiness);
   fChain->SetBranchAddress("fqnse", &fqnse, &b_fqnse);
   fChain->SetBranchAddress("fqitwnd", fqitwnd, &b_fqitwnd);
   fChain->SetBranchAddress("fqipeak", fqipeak, &b_fqipeak);
   fChain->SetBranchAddress("fqnhitpmt", fqnhitpmt, &b_fqnhitpmt);
   fChain->SetBranchAddress("fqtotq", fqtotq, &b_fqtotq);
   fChain->SetBranchAddress("fq0rtotmu", fq0rtotmu, &b_fq0rtotmu);
   fChain->SetBranchAddress("fq0rnll", fq0rnll, &b_fq0rnll);
   fChain->SetBranchAddress("fqn50", fqn50, &b_fqn50);
   fChain->SetBranchAddress("fqq50", fqq50, &b_fqq50);
   fChain->SetBranchAddress("fq1rpcflg", fq1rpcflg, &b_fq1rpcflg);
   fChain->SetBranchAddress("fq1rmom", fq1rmom, &b_fq1rmom);
   fChain->SetBranchAddress("fq1rt0", fq1rt0, &b_fq1rt0);
   fChain->SetBranchAddress("fq1rtotmu", fq1rtotmu, &b_fq1rtotmu);
   fChain->SetBranchAddress("fq1rnll", fq1rnll, &b_fq1rnll);
   fChain->SetBranchAddress("fq1rpos", fq1rpos, &b_fq1rpos);
   fChain->SetBranchAddress("fq1rdir", fq1rdir, &b_fq1rdir);
   fChain->SetBranchAddress("fq1rdconv", fq1rdconv, &b_fq1rdconv);
   fChain->SetBranchAddress("fq1reloss", fq1reloss, &b_fq1reloss);
   fChain->SetBranchAddress("fqpi0pcflg", fqpi0pcflg, &b_fqpi0pcflg);
   fChain->SetBranchAddress("fqpi0mom1", fqpi0mom1, &b_fqpi0mom1);
   fChain->SetBranchAddress("fqpi0mom2", fqpi0mom2, &b_fqpi0mom2);
   fChain->SetBranchAddress("fqpi0momtot", fqpi0momtot, &b_fqpi0momtot);
   fChain->SetBranchAddress("fqpi0dconv1", fqpi0dconv1, &b_fqpi0dconv1);
   fChain->SetBranchAddress("fqpi0dconv2", fqpi0dconv2, &b_fqpi0dconv2);
   fChain->SetBranchAddress("fqpi0t0", fqpi0t0, &b_fqpi0t0);
   fChain->SetBranchAddress("fqpi0totmu", fqpi0totmu, &b_fqpi0totmu);
   fChain->SetBranchAddress("fqpi0nll", fqpi0nll, &b_fqpi0nll);
   fChain->SetBranchAddress("fqpi0mass", fqpi0mass, &b_fqpi0mass);
   fChain->SetBranchAddress("fqpi0photangle", fqpi0photangle, &b_fqpi0photangle);
   fChain->SetBranchAddress("fqpi0pos", fqpi0pos, &b_fqpi0pos);
   fChain->SetBranchAddress("fqpi0dir1", fqpi0dir1, &b_fqpi0dir1);
   fChain->SetBranchAddress("fqpi0dir2", fqpi0dir2, &b_fqpi0dir2);
   fChain->SetBranchAddress("fqpi0dirtot", fqpi0dirtot, &b_fqpi0dirtot);
   fChain->SetBranchAddress("fqpmgpcflg", fqpmgpcflg, &b_fqpmgpcflg);
   fChain->SetBranchAddress("fqpmgmom1", fqpmgmom1, &b_fqpmgmom1);
   fChain->SetBranchAddress("fqpmgmom2", fqpmgmom2, &b_fqpmgmom2);
   fChain->SetBranchAddress("fqpmgt01", fqpmgt01, &b_fqpmgt01);
   fChain->SetBranchAddress("fqpmgt02", fqpmgt02, &b_fqpmgt02);
   fChain->SetBranchAddress("fqpmgtotmu", fqpmgtotmu, &b_fqpmgtotmu);
   fChain->SetBranchAddress("fqpmgnll", fqpmgnll, &b_fqpmgnll);
   fChain->SetBranchAddress("fqpmgpos1", fqpmgpos1, &b_fqpmgpos1);
   fChain->SetBranchAddress("fqpmgpos2", fqpmgpos2, &b_fqpmgpos2);
   fChain->SetBranchAddress("fqpmgdir1", fqpmgdir1, &b_fqpmgdir1);
   fChain->SetBranchAddress("fqpmgdir2", fqpmgdir2, &b_fqpmgdir2);
   fChain->SetBranchAddress("fqnmrfit", &fqnmrfit, &b_fqnmrfit);
   fChain->SetBranchAddress("fqmrifit", fqmrifit, &b_fqmrifit);
   fChain->SetBranchAddress("fqmrnring", fqmrnring, &b_fqmrnring);
   fChain->SetBranchAddress("fqmrpcflg", fqmrpcflg, &b_fqmrpcflg);
   fChain->SetBranchAddress("fqmrnll", fqmrnll, &b_fqmrnll);
   fChain->SetBranchAddress("fqmrtotmu", fqmrtotmu, &b_fqmrtotmu);
   fChain->SetBranchAddress("fqmrpid", fqmrpid, &b_fqmrpid);
   fChain->SetBranchAddress("fqmrmom", fqmrmom, &b_fqmrmom);
   fChain->SetBranchAddress("fqmrdconv", fqmrdconv, &b_fqmrdconv);
   fChain->SetBranchAddress("fqmreloss", fqmreloss, &b_fqmreloss);
   fChain->SetBranchAddress("fqmrt0", fqmrt0, &b_fqmrt0);
   fChain->SetBranchAddress("fqmrpos", fqmrpos, &b_fqmrpos);
   fChain->SetBranchAddress("fqmrdir", fqmrdir, &b_fqmrdir);
   fChain->SetBranchAddress("fqmsnfit", &fqmsnfit, &b_fqmsnfit);
   fChain->SetBranchAddress("fqmspcflg", fqmspcflg, &b_fqmspcflg);
   fChain->SetBranchAddress("fqmsnseg", fqmsnseg, &b_fqmsnseg);
   fChain->SetBranchAddress("fqmspid", fqmspid, &b_fqmspid);
   fChain->SetBranchAddress("fqmsifit", fqmsifit, &b_fqmsifit);
   fChain->SetBranchAddress("fqmsimer", fqmsimer, &b_fqmsimer);
   fChain->SetBranchAddress("fqmstotmu", fqmstotmu, &b_fqmstotmu);
   fChain->SetBranchAddress("fqmsnll", fqmsnll, &b_fqmsnll);
   fChain->SetBranchAddress("fqmsmom", fqmsmom, &b_fqmsmom);
   fChain->SetBranchAddress("fqmseloss", fqmseloss, &b_fqmseloss);
   fChain->SetBranchAddress("fqmst0", fqmst0, &b_fqmst0);
   fChain->SetBranchAddress("fqmspos", fqmspos, &b_fqmspos);
   fChain->SetBranchAddress("fqmsdir", fqmsdir, &b_fqmsdir);
   fChain->SetBranchAddress("fqtestn1r", &fqtestn1r, &b_fqtestn1r);
   fChain->SetBranchAddress("fqtest1rstage", fqtest1rstage, &b_fqtest1rstage);
   fChain->SetBranchAddress("fqtest1rse", fqtest1rse, &b_fqtest1rse);
   fChain->SetBranchAddress("fqtest1rpid", fqtest1rpid, &b_fqtest1rpid);
   fChain->SetBranchAddress("fqtest1rpcflg", fqtest1rpcflg, &b_fqtest1rpcflg);
   fChain->SetBranchAddress("fqtest1rmom", fqtest1rmom, &b_fqtest1rmom);
   fChain->SetBranchAddress("fqtest1rt0", fqtest1rt0, &b_fqtest1rt0);
   fChain->SetBranchAddress("fqtest1rtotmu", fqtest1rtotmu, &b_fqtest1rtotmu);
   fChain->SetBranchAddress("fqtest1rnll", fqtest1rnll, &b_fqtest1rnll);
   fChain->SetBranchAddress("fqtest1rpos", fqtest1rpos, &b_fqtest1rpos);
   fChain->SetBranchAddress("fqtest1rdir", fqtest1rdir, &b_fqtest1rdir);
   fChain->SetBranchAddress("fqtest1rdconv", fqtest1rdconv, &b_fqtest1rdconv);
   fChain->SetBranchAddress("fqtest1reloss", fqtest1reloss, &b_fqtest1reloss);
   fChain->SetBranchAddress("fqtestnpi0", &fqtestnpi0, &b_fqtestnpi0);
   fChain->SetBranchAddress("fqtestpi0stage", fqtestpi0stage, &b_fqtestpi0stage);
   fChain->SetBranchAddress("fqtestpi0pcflg", fqtestpi0pcflg, &b_fqtestpi0pcflg);
   fChain->SetBranchAddress("fqtestpi0mom1", fqtestpi0mom1, &b_fqtestpi0mom1);
   fChain->SetBranchAddress("fqtestpi0mom2", fqtestpi0mom2, &b_fqtestpi0mom2);
   fChain->SetBranchAddress("fqtestpi0momtot", fqtestpi0momtot, &b_fqtestpi0momtot);
   fChain->SetBranchAddress("fqtestpi0dconv1", fqtestpi0dconv1, &b_fqtestpi0dconv1);
   fChain->SetBranchAddress("fqtestpi0dconv2", fqtestpi0dconv2, &b_fqtestpi0dconv2);
   fChain->SetBranchAddress("fqtestpi0t0", fqtestpi0t0, &b_fqtestpi0t0);
   fChain->SetBranchAddress("fqtestpi0totmu", fqtestpi0totmu, &b_fqtestpi0totmu);
   fChain->SetBranchAddress("fqtestpi0nll", fqtestpi0nll, &b_fqtestpi0nll);
   fChain->SetBranchAddress("fqtestpi0mass", fqtestpi0mass, &b_fqtestpi0mass);
   fChain->SetBranchAddress("fqtestpi0photangle", fqtestpi0photangle, &b_fqtestpi0photangle);
   fChain->SetBranchAddress("fqtestpi0pos", fqtestpi0pos, &b_fqtestpi0pos);
   fChain->SetBranchAddress("fqtestpi0dir1", fqtestpi0dir1, &b_fqtestpi0dir1);
   fChain->SetBranchAddress("fqtestpi0dir2", fqtestpi0dir2, &b_fqtestpi0dir2);
   fChain->SetBranchAddress("fqtestpi0dirtot", fqtestpi0dirtot, &b_fqtestpi0dirtot);
   fChain->SetBranchAddress("prid_ang", &prid_ang, &b_prid_ang);
   fChain->SetBranchAddress("prid_len", &prid_len, &b_prid_len);
   fChain->SetBranchAddress("prid_pmom", &prid_pmom, &b_prid_pmom);
   fChain->SetBranchAddress("prid_mang", &prid_mang, &b_prid_mang);
   fChain->SetBranchAddress("prid_plike", &prid_plike, &b_prid_plike);
   fChain->SetBranchAddress("prid_mlike", &prid_mlike, &b_prid_mlike);
   fChain->SetBranchAddress("prid_nused", &prid_nused, &b_prid_nused);
   fChain->SetBranchAddress("dirtotepi", dirtotepi, &b_dirtotepi);
   fChain->SetBranchAddress("dirtotenpi", dirtotenpi, &b_dirtotenpi);
   fChain->SetBranchAddress("dirtotmue", dirtotmue, &b_dirtotmue);
   fChain->SetBranchAddress("dirsum", dirsum, &b_dirsum);
   fChain->SetBranchAddress("etot", &etot, &b_etot);
   fChain->SetBranchAddress("etotepi", &etotepi, &b_etotepi);
   fChain->SetBranchAddress("etotenpi", &etotenpi, &b_etotenpi);
   fChain->SetBranchAddress("etotmue", &etotmue, &b_etotmue);
   fChain->SetBranchAddress("oscweight", oscweight, &b_oscweight);
   fChain->SetBranchAddress("oscwgt", &oscwgt, &b_oscwgt);
   fChain->SetBranchAddress("live", &live, &b_live);
   fChain->SetBranchAddress("sacth", &sacth, &b_sacth);
   fChain->SetBranchAddress("sactg", &sactg, &b_sactg);
   fChain->SetBranchAddress("sacth1d", &sacth1d, &b_sacth1d);
   fChain->SetBranchAddress("ent_pos", ent_pos, &b_ent_pos);
   fChain->SetBranchAddress("ent_dir", ent_dir, &b_ent_dir);
   fChain->SetBranchAddress("length", &length, &b_length);
   fChain->SetBranchAddress("tr_mom1", &tr_mom1, &b_tr_mom1);
   fChain->SetBranchAddress("A_ent_mom", &A_ent_mom, &b_A_ent_mom);
   fChain->SetBranchAddress("A_ent_pos", A_ent_pos, &b_A_ent_pos);
   fChain->SetBranchAddress("A_ent_dir", A_ent_dir, &b_A_ent_dir);
   fChain->SetBranchAddress("A_ext_mom", &A_ext_mom, &b_A_ext_mom);
   fChain->SetBranchAddress("A_ext_pos", A_ext_pos, &b_A_ext_pos);
   fChain->SetBranchAddress("A_ext_dir", A_ext_dir, &b_A_ext_dir);
   fChain->SetBranchAddress("Fit_pos", Fit_pos, &b_Fit_pos);
   fChain->SetBranchAddress("Fit_dir", Fit_dir, &b_Fit_dir);
   fChain->SetBranchAddress("Fit_len", &Fit_len, &b_Fit_len);
   fChain->SetBranchAddress("Fit_mom", &Fit_mom, &b_Fit_mom);
   fChain->SetBranchAddress("Fit_pid", &Fit_pid, &b_Fit_pid);
   fChain->SetBranchAddress("Um_ehit8m", &Um_ehit8m, &b_Um_ehit8m);
   fChain->SetBranchAddress("Um_ohit8m", &Um_ohit8m, &b_Um_ohit8m);
   fChain->SetBranchAddress("Um_qent", &Um_qent, &b_Um_qent);
   fChain->SetBranchAddress("Sh_chi1p", &Sh_chi1p, &b_Sh_chi1p);
   fChain->SetBranchAddress("Sh_delta", &Sh_delta, &b_Sh_delta);
   fChain->SetBranchAddress("Sh_mean", &Sh_mean, &b_Sh_mean);
   fChain->SetBranchAddress("Sh_meanq", &Sh_meanq, &b_Sh_meanq);
   fChain->SetBranchAddress("Sh_id", &Sh_id, &b_Sh_id);
   fChain->SetBranchAddress("Tr_stop", Tr_stop, &b_Tr_stop);
   fChain->SetBranchAddress("Tr_mom", &Tr_mom, &b_Tr_mom);
   fChain->SetBranchAddress("Tr_len", &Tr_len, &b_Tr_len);
   fChain->SetBranchAddress("Tr_len1", &Tr_len1, &b_Tr_len1);
   fChain->SetBranchAddress("Pid_flg", &Pid_flg, &b_Pid_flg);
   fChain->SetBranchAddress("Crs1", &Crs1, &b_Crs1);
   fChain->SetBranchAddress("Crs2", &Crs2, &b_Crs2);
   fChain->SetBranchAddress("iclass", &iclass, &b_iclass);
   fChain->SetBranchAddress("mu_class", &mu_class, &b_mu_class);
   fChain->SetBranchAddress("mu_dec", &mu_dec, &b_mu_dec);
   fChain->SetBranchAddress("mu_dir", mu_dir, &b_mu_dir);
   fChain->SetBranchAddress("mu_pos", mu_pos, &b_mu_pos);
   fChain->SetBranchAddress("mu_good", &mu_good, &b_mu_good);
   fChain->SetBranchAddress("history", &history, &b_history);
   fChain->SetBranchAddress("Pdst", &Pdst, &b_Pdst);
   fChain->SetBranchAddress("idoff", &idoff, &b_idoff);
   fChain->SetBranchAddress("anthit", &anthit, &b_anthit);
   fChain->SetBranchAddress("idseq", &idseq, &b_idseq);
   fChain->SetBranchAddress("tstfrac", &tstfrac, &b_tstfrac);
   fChain->SetBranchAddress("judge", &judge, &b_judge);
   fChain->SetBranchAddress("Upcrs1", &Upcrs1, &b_Upcrs1);
   fChain->SetBranchAddress("Upcrs2", &Upcrs2, &b_Upcrs2);
   fChain->SetBranchAddress("lst", &lst, &b_lst);
   fChain->SetBranchAddress("jd", &jd, &b_jd);
   fChain->SetBranchAddress("fjd", &fjd, &b_fjd);
   fChain->SetBranchAddress("alt", &alt, &b_alt);
   fChain->SetBranchAddress("azi", &azi, &b_azi);
   fChain->SetBranchAddress("ra", &ra, &b_ra);
   fChain->SetBranchAddress("dec", &dec, &b_dec);
   fChain->SetBranchAddress("glat", &glat, &b_glat);
   fChain->SetBranchAddress("glong", &glong, &b_glong);
   Notify();
}

Bool_t NewVars::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NewVars::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NewVars::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NewVars_cxx
