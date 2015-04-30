#include "StringArray.h"
#include "Error.h"
#include <stdio.h>
#include "math.h"
#include <iostream>
#include "StringHash.h"
#include "MathStats.h"
#include "MixChidist.h"
#include "MathSVD.h"
#include "Meta.h"
#define MATHLIB_STANDALONE 
#include <Rmath.h>
#include "Calculate_mvt_pvalue.h"
#include "My_mvt.h"
#include "InputFile.h"
#include "SummaryFileReader.h"
#include "QuickIndex.h"
#include "cquad.h"
#include <Eigen/Eigenvalues>

String Meta::studyName = "";
String Meta::prefix = ""; 
bool Meta::correctGC = false;
bool Meta::Burden = false;
bool Meta::MB = false;
bool Meta::SKAT = false;
bool Meta::SKATO = false;
bool Meta::VTa = false;
bool Meta::VTp = false;
bool Meta::outvcf =false; 
bool Meta::fullResult =false; 
bool Meta::report = false; 
double  Meta::report_pvalue_cutoff = 1e-06; 
bool Meta::founderAF = false; 
double Meta::HWE = 1e-05;
double Meta::CALLRATE = 0.95;
double Meta::MAF_cutoff = 0.05;
int Meta::marker_col = 2;
int Meta::cov_col = 3;
bool Meta::tabix = false;
String Meta::cond = "";

Meta::Meta(){}
Meta::~Meta()
{}

//This function read all study names from a file 
//and save the names in the StringArray files
void Meta::Prepare()
{
   if(prefix != "")
   {
      if(prefix.Last()=='.' || prefix.Last()=='/')
	 pdf_filename = prefix + "meta.plots.pdf";
      else
	 pdf_filename = prefix + ".meta.plots.pdf";
   }
   else
      pdf_filename = "meta.plots.pdf";
   pdf.OpenFile(pdf_filename);


   if(studyName!="")
   {
      IFILE inFile = ifopen(studyName,"r");
      if(inFile==NULL)
	 error("FATAL ERROR! Please check file name for --studyName  %s\n",studyName.c_str());
      String buffer;
      while (!ifeof(inFile))
      {
	 buffer.ReadLine(inFile);
	 if(buffer.FindChar('#')!=-1)
	    continue;
	 studies.AddTokens(buffer, "\n");
      }
      ifclose(inFile);
   }
   else
      error("FATAL ERROR! --studyName can not be empty! Usage: --studyName yourStudynames\n");

   //if conditional analysis says yes
   if(cond!="")
   {
      IFILE condFile = ifopen(cond,"r");
      if(condFile==NULL)
	 error("Can not open file %s.\n",cond.c_str());
      String buffer;
      while (!ifeof(condFile))
      {
	 buffer.ReadLine(condFile);
	 if(buffer.FindChar('#')!=-1 || buffer.FindChar(':')==-1)
	    continue;
	 StringArray tmpMarker;
	 tmpMarker.AddTokens(buffer, " \t");
	 for(int i=0;i<tmpMarker.Length();i++)
	 {
	    commonVar.Push(tmpMarker[i]);
	    StringArray tmp;
	    tmp.AddTokens(tmpMarker[i],":");
	    common_chr.Push(tmp[0]);
	    common_pos.Push(tmp[1]);
	    common_ref.Push(tmp[2]);
	    common_alt.Push(tmp[3]);
	    conditionVar.SetInteger(tmp[0]+":"+ tmp[1],commonVar.Length()-1);
	 }
      }
      ifclose(condFile);
      if(common_pos.Length()==0)
	 error("No variant to be conditioned on. Condition file %s might be in the wrong format.\n",cond.c_str());

      int numStudy = studies.Length();
      commonVar_study = new IntArray [numStudy];
      commonVar_effSize = new Vector [numStudy];
      commonVar_U = new Vector [numStudy];
      commonVar_V = new Vector [numStudy];
      commonVar_markers_in_window = new IntArray * [numStudy];
      commonVar_marker_cov = new Vector * [numStudy];
      XX_inv = new Matrix [numStudy];
      cond_status = new bool [numStudy];
      commonVar_betaHat = new Vector [numStudy];

      for(int s=0;s<studies.Length();s++)
      {
	 SummaryFileReader statReader,covReader;

	 String filename = studies[s] + ".singlevar.score.txt.gz";
	 String covFilename = studies[s] + ".singlevar.cov.txt.gz";
	 int adjust =0;
	 if(!statReader.ReadTabix(filename))
	 {
	    adjust =1;
	    filename = studies[s] + ".MetaScore.assoc.gz";
	    if(!statReader.ReadTabix(filename))
	    {
	       error("Cannot open file: %s.singlevar.score.txt.gz!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip yourfile.singlevar.score.txt; tabix -c \"#\" -s 1 -b 2 -e 2 yourprefix.singlevar.score.txt.gz\n",studies[s].c_str());
	    }
	 }

	 IFILE dup = ifopen(filename,"r");
	 buffer.ReadLine(dup);
	 if(buffer.Find("RareMetalWorker")==-1)
	    adjust =1;
	 ifclose(dup);

	 dup=ifopen(filename,"r"); 

	 while (!ifeof(dup))
	 {
	    buffer.ReadLine(dup);
	    //get sample size
	    if(buffer.Find("##AnalyzedSamples")!=-1)
	    {
	       StringArray tokens;
	       tokens.AddTokens(buffer,"=");
	       SampleSize.Push(tokens[1].AsInteger());
	       break;
	    }
	 }
	 ifclose(dup);

	 if(!covReader.ReadTabix(covFilename))
	 {
	    covFilename = studies[s] + ".MetaCov.assoc.gz";
	    if(!covReader.ReadTabix(covFilename))
	    {
	       error("Cannot open file: %s.singlevar.cov.txt.gz!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip yourfile.singlevar.score.txt; tabix -c \"#\" -s 1 -b 2 -e 2 yourprefix.singlevar.score.txt.gz\n",studies[s].c_str());
	    }
	 }

	 for(int i=0;i<commonVar.Length();i++)
	 {
	    //if this variant is genotyped in this study
	    if(!covReader.ReadRecord(common_chr[i],common_pos[i]))
	       continue;
	    if(!statReader.ReadRecord(common_chr[i],common_pos[i]))
	       continue;
	    StringArray record;
	    record.AddTokens(statReader.buffer,"\t");

	    if((record[2]==common_ref[i] && record[3]==common_alt[i])
		  || (record[3]==common_ref[i] && record[2]==common_alt[i]))
	    {
	       double v = record[14-adjust].AsDouble();
	       if(v>0)
	       {
		  commonVar_study[s].Push(i);
		  commonVar_V[s].Push(v*v);
		  commonVar_effSize[s].Push(record[15-adjust].AsDouble());
		  commonVar_U[s].Push(record[13-adjust].AsDouble());
	       }
	    }
	 }
	 int dim = commonVar_study[s].Length();
	 if(dim!=0)
	    cond_status[s] =true;
	 else 
	 {
	    cond_status[s] = false;
	    continue;
	 }
	 commonVar_markers_in_window[s] = new IntArray [dim];
	 commonVar_marker_cov[s] = new Vector [dim];

	 covReader.ReadTabix(covFilename);
	 for(int i=0;i<dim;i++)
	 {
	    int idx = commonVar_study[s][i];
	    if(!covReader.ReadRecord(common_chr[idx],common_pos[idx]))
	       continue;
	    StringArray tmp_markers;
	    tmp_markers.AddTokens(covReader.marker_nearby,",");
	    for(int j=0;j<tmp_markers.Length();j++)
	       commonVar_markers_in_window[s][i].Push(tmp_markers[j].AsInteger());
	    tmp_markers.Clear();
	    tmp_markers.AddTokens(covReader.marker_cov,",");
	    for(int j=0;j<tmp_markers.Length();j++)
	       commonVar_marker_cov[s][i].Push(tmp_markers[j].AsDouble());
	 }
	 XX_inv[s].Dimension(commonVar_study[s].Length(),commonVar_study[s].Length(),0.0);
	 CalculateXXCov(s,XX_inv[s]);
	 //calculate beta hat for common variants
	 commonVar_betaHat[s].Dimension(dim);
	 for(int i=0;i<commonVar_betaHat[s].dim;i++)
	    commonVar_betaHat[s][i] = XX_inv[s][i].InnerProduct(commonVar_U[s]);
	 //printf("xx_inv %g %g %g,U %g %g\n",XX_inv[s][0][0],XX_inv[s][0][1],XX_inv[s][1][1],commonVar_U[0][0],commonVar_U[0][1]);
	 //printf("xx_inv %g,commonvar U %g commonVar_betahat %g\n",XX_inv[s][0][0],commonVar_U[s][0],commonVar_betaHat[s][0]);
      }

      bool status = false;
      for(int i=0;i<numStudy;i++)
      {
	 if(cond_status[i]) 
	 {
	    status=true;
	    break;
	 }
      }
      if(!status) 
      {
	 cond="";
	 printf("\nWarning: None of the variants to be conditioned on are found in any of the studies. Conditional analysis options will be ignored.\n\n");
      }
   }
}

//this function will read through summary statistics of each study
//and pool the information. 
//At the end, single variant meta-analysis will be completed.
void Meta::PoolSummaryStat(GroupFromAnnotation & group)
{
   //usefulSize and usefulAC have the pooled N and AC information 
   StringIntHash usefulSize;
   StringDoubleHash usefulAC; 
   //refalt use study:chr:pos as key and ref:alt as value
   StringHash refalt;
   StringIntHash SNP_DirectionByStudy,SNP_DirectionByStudyNoAllele;
   StringArray directions;

   Vector GCbyStudy;
   total_N=0;
   int flipCount=0,skip_count=0;

   for(int study=0;study<studies.Length();study++)
   {
      Vector chisq_study_i;
      //#CHROM  POS     REF     ALT     N_INFORMATIVE   FOUNDER_AF      ALL_AF  INFORMATIVE_ALT_AC      CALL_RATE       HWE_PVALUE      N_REF   N_HET   N_ALT   U_STAT  SQRT_V_STAT     ALT_EFFSIZE     PVALUE

      printf("Pooling summary statistics from study %s ...\n",studies[study].c_str());

      //read in summary statistics. 
      //maf and summary stat are saved. SNPpos in the file are hashed.
      String filename = studies[study] + ".singlevar.score.txt.gz";

      SummaryFileReader covReader;
      if(cond!="" && cond_status[study])
      {
	 String covFilename = studies[study] + ".singlevar.cov.txt.gz";
	 covReader.ReadTabix(covFilename);
      }

      IFILE file;
      file = ifopen(filename,"r");

      int adjust=0;

      if(file == NULL)
      {
	 filename = studies[study] + ".MetaScore.assoc.gz";
	 file = ifopen(filename,"r");
	 adjust=1;
      }

      if(file == NULL)
	 error("Cannot open file: %s.singlevar.score.txt.gz!\nInput file has to be bgzipped and tabix indexed using the following command:\n bgzip yourfile.singlevar.score.txt; tabix -c \"#\" -s 1 -b 2 -e 2 yourprefix.singlevar.score.txt.gz\n",studies[study].c_str());

      String buffer;
      StringArray tokens;
      int current_N =0;

      buffer.ReadLine(file);
      if(buffer.Find("RareMetalWorker")==-1)
	 adjust =1;
      ifclose(file);
      FormatAdjust.Push(adjust);

      file = ifopen(filename,"r");
      while (!ifeof(file))
      {
	 buffer.ReadLine(file);
	 //get sample size
	 if(buffer.Find("##AnalyzedSamples")!=-1)
	 {
	    StringArray tokens;
	    tokens.AddTokens(buffer,"=");
	    SampleSize.Push(tokens[1].AsInteger());
	    continue;
	 }
	 else if(buffer.FindChar('#') !=-1 || buffer.Find("CHROM")!=-1)
	    continue;
	 tokens.Clear();
	 tokens.AddTokens(buffer, SEPARATORS);

	 if(tokens[0].Find("chr")!=-1)
	    tokens[0] = tokens[0].SubStr(3);

	 //push the chisq statistics
	 if(tokens[14-adjust].AsDouble()>0)
	 {
	    double u = tokens[13-adjust].AsDouble();
	    double v = tokens[14-adjust].AsDouble();
	    chisq_study_i.Push(u*u/(v*v));
	 }

	 //exclude the variants that have one allele missing or are monomorphic
	 if(tokens[2]=="." || tokens[3]=="." || tokens[2]==tokens[3] || tokens[5].AsDouble()==1.0 || tokens[5].AsDouble()==0.0) 
	    continue;

	 String tmp_SNP = tokens[0]+":"+tokens[1]+":"+tokens[2]+":"+tokens[3];
	 String tmp_SNP_flip = tokens[0]+":"+tokens[1]+":"+tokens[3]+":"+tokens[2];
	 String refalt_current = study + ":" + tokens[0] + ":" + tokens[1];

	 //check if this SNP (allowing allele flips) has been hashed before
	 int index1 = usefulAC.Find(tmp_SNP);
	 int index2 = usefulAC.Find(tmp_SNP_flip);

	 //QC: filtering variants
	 if(tokens[8-adjust].AsDouble()<CALLRATE || tokens[9-adjust].AsDouble()<HWE)
	 {
	    int directions_idx;
	    if(index1 !=-1)
	    {
	       directions_idx = SNP_DirectionByStudy.Integer(tmp_SNP);
	       if(directions[directions_idx].Length()<study)
	       {
		  int len = directions[directions_idx].Length();
		  for(int s=0;s<study-len;s++)
		     directions[directions_idx] += "?";
	       }
	       directions[directions_idx] += "!";
	    }
	    else if(index2 != -1)
	    {
	       directions_idx = SNP_DirectionByStudy.Integer(tmp_SNP_flip);
	       if(directions[directions_idx].Length()<study)
	       {
		  int len = directions[directions_idx].Length();
		  for(int s=0;s<study-len;s++)
		     directions[directions_idx] += "?";
	       }
	       directions[directions_idx] += "!";
	    }
	    else if(tokens[2]!=tokens[3])
	    {
	       String * refalt_saved = (String *) refalt.Object(tokens[0]+":"+tokens[1]);
	       if(refalt_saved != NULL)
	       {
		  //if this variant at this chr:pos has been hashed but with different ref/alt alleles
		  //then skip this variant
		  skip_count++;
		  String skip_SNP;
		  skip_SNP = study;
		  skip_SNP += ":";
		  skip_SNP +=tokens[0];
		  skip_SNP += ":";
		  skip_SNP +=tokens[1];

		  SNPexclude.SetInteger(skip_SNP,skip_count);
		  //printf("this SNP is being skipped %d %s %s %d \n",study,tokens[0].c_str(),tokens[1].c_str(),skip_count);
		  directions_idx = SNP_DirectionByStudyNoAllele.Integer(tokens[0]+":"+tokens[1]);
		  if(directions[directions_idx].Length()<study)
		  {
		     int len = directions[directions_idx].Length();
		     for(int s=0;s<study-len;s++)
			directions[directions_idx] += "?";
		  }
		  directions[directions_idx] += "!";
		  continue;
	       }

	       String direction = "";
	       if(study>0)
	       {
		  for(int s=0;s<study;s++)
		     direction += "?";
	       }
	       direction += "!";
	       directions.Push(direction);
	       SNP_DirectionByStudy.SetInteger(tmp_SNP,directions.Length()-1);
	       SNP_DirectionByStudyNoAllele.SetInteger(tmp_SNP,directions.Length()-1);
	    }
	    continue;
	 }

	 int current_AC;
	 if(tokens[5]=="NA")
	 {
	    current_AC = 0;
	    current_N = SampleSize[study];
	 }
	 else
	 {
	    current_AC = 2*tokens[12-adjust].AsInteger() + tokens[11-adjust].AsInteger();
	    current_N = tokens[10-adjust].AsInteger()+tokens[11-adjust].AsInteger()+tokens[12-adjust].AsInteger();
	 }

	 if(index1 !=-1 || index2 !=-1)
	 {
	    if(index2!=-1)
	    {
	       flipCount++;
	       flipSNP.SetInteger(study+":"+tokens[0]+":"+tokens[1],flipCount);
	    }
	    //update usefulAC
	    //handling flipped alleles
	    if(index2!=-1)	       
	    {
	       double AC = usefulAC.Double(index2);
	       if(founderAF)
	       {
		  AC += 2.0*SampleSize[study]*(1.0-tokens[5].AsDouble());
		  usefulAC.SetDouble(index2,AC);
	       }
	       else
	       {
		  AC += current_N*2.0-current_AC;  
		  usefulAC.SetDouble(index2,AC);
	       }
	    }
	    else 
	    {
	       double AC = usefulAC.Double(index1);
	       if(founderAF)
	       {
		  AC += 2.0*SampleSize[study]*tokens[5].AsDouble();
		  usefulAC.SetDouble(index1,AC);
	       }
	       else
	       {
		  AC += current_AC;
		  usefulAC.SetDouble(index1,AC);
	       }
	    }

	    //update usefulSize
	    int N_idx;
	    if(index1!=-1)
	       N_idx= usefulSize.Find(tmp_SNP);
	    else
	       N_idx= usefulSize.Find(tmp_SNP_flip);

	    int N = usefulSize.GetCount(N_idx); 

	    if(founderAF)
	    {
	       N+=SampleSize[study];
	    }
	    else
	    {
	       N += current_N;
	    }

	    usefulSize.SetInteger(N_idx,N);

	    //update SNPstat
	    int stat_idx;
	    if(index1!=-1)
	    {
	       stat_idx = SNPstat.Find(tmp_SNP);
	       double stat = SNPstat.Double(stat_idx);
	       stat += tokens[13-adjust].AsDouble();
	       SNPstat.SetDouble(stat_idx,stat);
	    }
	    else
	    {
	       stat_idx = SNPstat.Find(tmp_SNP_flip);
	       double stat = SNPstat.Double(stat_idx);
	       stat += -1.0*tokens[13-adjust].AsDouble();
	       SNPstat.SetDouble(stat_idx,stat);
	    }

	    //update SNP_Vstat
	    int Vstat_idx;
	    if(index1!=-1)
	       Vstat_idx = SNP_Vstat.Find(tmp_SNP);
	    else
	       Vstat_idx = SNP_Vstat.Find(tmp_SNP_flip);

	    double Vstat = SNP_Vstat.Double(Vstat_idx);
	    double vstat_tmp = tokens[14-adjust].AsDouble();
	    vstat_tmp = vstat_tmp*vstat_tmp;
	    Vstat += vstat_tmp;
	    SNP_Vstat.SetDouble(Vstat_idx,Vstat);

	    //update SNP_cond_stat and SNP_cond_V
	    if(cond!="" && cond_status[study] && tokens[14-adjust].AsDouble()>0.0)
	    {
	       if(conditionVar.Integer(tokens[0]+":"+tokens[1])==-1)
	       {
		  Vector GX;
		  CalculateGenotypeCov(covReader,tokens[0],tokens[1].AsInteger(),study,GX);      
		  double cond_u = tokens[13-adjust].AsDouble() - GX.InnerProduct(commonVar_betaHat[study]);
		  int idx;
		  if(index1!=-1)
		  {
		     idx = SNPstat_cond.Find(tmp_SNP);
		     double stat = SNPstat_cond.Double(idx);
		     stat += cond_u;
		     SNPstat_cond.SetDouble(idx,stat);
		  }
		  else
		  {
		     idx = SNPstat_cond.Find(tmp_SNP_flip);
		     double stat = SNPstat_cond.Double(idx);
		     stat += -1.0*cond_u;
		     SNPstat_cond.SetDouble(idx,stat);
		  }

		  Vector tmp;
		  for(int i=0;i<GX.dim;i++)
		     tmp.Push(GX.InnerProduct(XX_inv[study][i]));
		  double v = tokens[14-adjust].AsDouble();
		  double cond_v_part = tmp.InnerProduct(GX);
		  double cond_V = v*v-cond_v_part;
		  if(index1!=-1)
		     idx = SNP_Vstat_cond.Find(tmp_SNP);
		  else
		     idx = SNP_Vstat_cond.Find(tmp_SNP_flip);
		  double cond_Vstat = SNP_Vstat_cond.Double(idx);
		  cond_Vstat += cond_V;
		  SNP_Vstat_cond.SetDouble(idx,cond_Vstat);
	       }
	    }
	    //update SNP_DirectionByStudy
	    int directions_idx;
	    String direction = "+";
	    if( tokens[15-adjust]=="-nan" || tokens[15-adjust]=="NA" || tokens[15-adjust]=="na" || tokens[15-adjust]=="nan" || tokens[15-adjust]=="-")
	    {
	       direction = "?";
	    }
	    if(index1!=-1)
	    {
	       if(tokens[15-adjust].AsDouble()<0)
		  direction = "-";
	       directions_idx = SNP_DirectionByStudy.Integer(tmp_SNP);	
	    }
	    else
	    {
	       if(tokens[15-adjust].AsDouble()>0)
		  direction = "-";
	       directions_idx = SNP_DirectionByStudy.Integer(tmp_SNP_flip);	
	    }
	    if(directions[directions_idx].Length()<study)
	    {
	       int len = directions[directions_idx].Length();
	       for(int s=0;s<study-len;s++)
		  directions[directions_idx] += "?";
	    }
	    directions[directions_idx] += direction;
	 }
	 else
	 {
	    String * refalt_saved = (String *) refalt.Object(tokens[0]+":"+tokens[1]);
	    if(refalt_saved != NULL)
	    {
	       //if this variant at this chr:pos has been hashed but with different ref/alt alleles
	       //then skip this variant
	       skip_count++;
	       String skip_SNP;
	       skip_SNP = study;
	       skip_SNP += ":";
	       skip_SNP +=tokens[0];
	       skip_SNP += ":";
	       skip_SNP +=tokens[1];

	       SNPexclude.SetInteger(skip_SNP,skip_count);
	       //printf("this SNP is being skipped %d %s %s %d \n",study,tokens[0].c_str(),tokens[1].c_str(),skip_count);
	       int directions_idx = SNP_DirectionByStudyNoAllele.Integer(tokens[0]+":"+tokens[1]);
	       if(directions[directions_idx].Length()<study)
	       {
		  int len = directions[directions_idx].Length();
		  for(int s=0;s<study-len;s++)
		     directions[directions_idx] += "?";
	       }
	       directions[directions_idx] += "!";
	       continue;
	    }

	    if(founderAF)
	    {
	       usefulAC.SetDouble(tmp_SNP,2.0*SampleSize[study]*tokens[5].AsDouble());
	       usefulSize.SetInteger(tmp_SNP,SampleSize[study]);
	    }
	    else
	    {
	       usefulAC.SetDouble(tmp_SNP,current_AC);
	       usefulSize.SetInteger(tmp_SNP,current_N);
	    }
	    SNPstat.SetDouble(tmp_SNP,tokens[13-adjust].AsDouble());
	    double vstat_tmp = tokens[14-adjust].AsDouble();
	    vstat_tmp = vstat_tmp*vstat_tmp;
	    SNP_Vstat.SetDouble(tmp_SNP,vstat_tmp);

	    if(cond!="" && cond_status[study])
	    {
	       Vector GX;
	       CalculateGenotypeCov(covReader,tokens[0],tokens[1].AsInteger(),study,GX);
	       double cond_u = tokens[13-adjust].AsDouble() - GX.InnerProduct(commonVar_betaHat[study]);
	       Vector tmp;
	       for(int i=0;i<GX.dim;i++)
		  tmp.Push(GX.InnerProduct(XX_inv[study][i]));
	       double cond_V = vstat_tmp - tmp.InnerProduct(GX);
	       //printf("N %d, SNP %s, U %g, beta %g, GX %g,cond_u %g, cond_v %g\n",tokens[4-adjust].AsInteger(),tmp_SNP.c_str(),tokens[13-adjust].AsDouble(),commonVar_betaHat[study][0],GX[0],cond_u, cond_V);
	       SNPstat_cond.SetDouble(tmp_SNP,cond_u);
	       SNP_Vstat_cond.SetDouble(tmp_SNP,cond_V);
	    }

	    refalt.SetObject(tokens[0]+":"+tokens[1], &refalt_current);

	    String direction = "";
	    if(study!=0)
	    {
	       for(int s=0;s<study;s++)
		  direction += "?";
	    }

	    if(tokens[15-adjust]=="NA" || tokens[15-adjust]=="na" || tokens[15-adjust]=="nan" || tokens[15-adjust]=="-nan"|| tokens[15-adjust]=="-")
	       direction += "?";
	    else if(tokens[15-adjust].AsDouble()<0)
	       direction += "-";
	    else if(tokens[15-adjust].AsDouble()>0)
	       direction += "+";

	    directions.Push(direction);
	    SNP_DirectionByStudy.SetInteger(tmp_SNP,directions.Length()-1);
	    SNP_DirectionByStudyNoAllele.SetInteger(tokens[0]+":"+tokens[1],directions.Length()-1);
	 }
      }
      ifclose(file);

      total_N += SampleSize[study];
      //calculate GC
      chisq_study_i.Sort();
      GCbyStudy.Push(chisq_study_i[0.5]/0.456);
   }

   //calculate pooled allele frequencies
   StringArray chr_AC,unique_chr,SNPname_AC;
   IntArray pos_AC;
   //get the unique chromosomes
   for(int i=0;i<usefulAC.Capacity();i++)
   {
      if(!usefulAC.SlotInUse(i))
      {
	 continue;
      }
      String SNPname = usefulAC[i];
      StringArray tmp;
      tmp.AddTokens(SNPname,":");
      chr_AC.Push(tmp[0]);
      pos_AC.Push(tmp[1].AsInteger());
      SNPname_AC.Push(SNPname);
      if(unique_chr.Find(tmp[0])==-1)
	 unique_chr.Push(tmp[0]);
   }
   QuickIndex chr_AC_idx(chr_AC);
   unique_chr.Sort();
   StringArray chr_cp,character_chr;
   for(int i=0;i<unique_chr.Length();i++)
   {
      if(unique_chr[i].AsInteger()<=22 && unique_chr[i].AsInteger()>=1)
	 chr_cp.Push(unique_chr[i]);
      else
	 character_chr.Push(unique_chr[i]);
   }
   for(int i=0;i<character_chr.Length();i++)
      chr_cp.Push(character_chr[i]);
   unique_chr = chr_cp; //now unique_chr are sorted as 1,2,...,22,X,Y,M...
   chr_cp.Clear();
   character_chr.Clear();
   for(int i=0;i<unique_chr.Length();i++)
   {
      IntArray pos_i;
      StringArray SNPname_i;
      for(int j=0;j<chr_AC.Length();j++)
      {
	 if(chr_AC[chr_AC_idx[j]] == unique_chr[i])
	 {
	    pos_i.Push(pos_AC[chr_AC_idx[j]]);
	    SNPname_i.Push(SNPname_AC[chr_AC_idx[j]]);
	 }
      }
      QuickIndex pos_i_idx(pos_i);
      for(int j=0;j<pos_i.Length();j++)
      {
	 double AC = usefulAC.Double(SNPname_i[pos_i_idx[j]]);
	 int N = usefulSize.Integer(SNPname_i[pos_i_idx[j]]);
	 double maf;
	 if(founderAF)
	 {
	    maf = AC/(2.0*total_N);
	 }
	 else
	 {
	    maf = AC/(2.0*N);
	 }
	 //printf("total_N is: %d\n N is: %d\n",total_N,N);
	 if(maf>0.5)
	    maf = 1.0-maf;

	 SNPmaf_maf.Push(maf);
	 SNPmaf_name.Push(SNPname_i[pos_i_idx[j]]);
	 SNPmaf.SetInteger(SNPname_i[pos_i_idx[j]],SNPmaf_maf.Length()-1);
      }
   }

   //calculate single variant meta results
   printf("\nPerforming Single variant meta analysis ...\n");
   //calculate final results here
   IFILE output;
   String filename;
   if(prefix =="")
      filename = "meta.singlevar.results";
   else if(prefix.Last()=='.' || prefix.Last()=='/')
      filename = prefix +  "meta.singlevar.results";
   else
      filename = prefix + ".meta.singlevar.results";

   output=ifopen(filename,"w",InputFile::UNCOMPRESSED);

   ifprintf(output,"##Method=SinglevarScore\n");
   ifprintf(output,"##STUDY_NUM=%d\n",studies.Length());
   ifprintf(output,"##TotalSampleSize=%d\n",total_N);
   if(cond=="")
      ifprintf(output,"#CHROM\tPOS\tREF\tALT\tPOOLED_ALT_AF\tDIRECTION_BY_STUDY\tEFFECT_SIZE\tPVALUE\n");
   else
      ifprintf(output,"#CHROM\tPOS\tREF\tALT\tPOOLED_ALT_AF\tDIRECTION_BY_STUDY\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");

   String vcf_filename;
   IFILE vcfout;
   if(outvcf)
   {
      if(prefix=="")
	 vcf_filename = "pooled.variants.vcf";
      else if(prefix.Last()=='.' || prefix.Last()=='/')
	 vcf_filename = prefix + "pooled.variants.vcf";
      else 
	 vcf_filename = prefix + ".pooled.variants.vcf";

      vcfout = ifopen(vcf_filename,"w",InputFile::UNCOMPRESSED);
      ifprintf(vcfout,"#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
   }

   Vector pvalueAll,pvalue1,pvalue5,pvalueAll_cond,pvalue1_cond,pvalue5_cond;
   StringArray chr_plot,geneLabel;
   Vector pos_plot,chisq_before_GC;
   IntArray pvalue1_idx,pvalue5_idx;

   //for annotation purpose
   String target_chr;
   int target_pos,target;
   double target_pvalue;
   //Sort variants by chr and pos
   for(int i=0;i<SNPmaf_maf.Length();i++)
   {
      String SNPname = SNPmaf_name[i];
      double U = SNPstat.Double(SNPname);
      double V = SNP_Vstat.Double(SNPname);
      double maf = SNPmaf_maf[i];
      int direction_idx = SNP_DirectionByStudy.Integer(SNPname);
      if(directions[direction_idx].Length()<studies.Length())
      {
	 int len = directions[direction_idx].Length();
	 for(int s=0;s<studies.Length()-len;s++)
	    directions[direction_idx] += "?";
      }

      String direction = directions[direction_idx];
      //printf("direction is %s,SNPname is %s\n",directions[direction_idx].c_str(), SNPname.c_str());

      if(maf>0.0 && maf<1.0 && V!=0.0)
      {
	 double chisq = U*U/V;
	 chisq_before_GC.Push(chisq);
	 double pvalue = pchisq(chisq,1,0,0);
	 double effSize = U/V;
	 double cond_pvalue = _NAN_,cond_effSize = _NAN_;
	 bool disect=false;
	 while(pvalue==0.0)
	 {
	    disect=true;
	    chisq *= 0.999;
	    pvalue = pchisq(chisq,1,0,0);
	 }
	 StringArray tokens;
	 tokens.AddTokens(SNPname,":");

	 singleVarPvalue.SetDouble(SNPname,pvalue);
	 singleVarEff.SetDouble(SNPname,effSize);

	 pvalueAll.Push(pvalue);
	 if(maf<0.01)
	 { 
	    pvalue1.Push(pvalue);
	    pvalue1_idx.Push(pvalueAll.Length()-1);
	 }
	 if(maf<0.05)
	 {
	    pvalue5.Push(pvalue);
	    pvalue5_idx.Push(pvalueAll.Length()-1);
	 }
	 chr_plot.Push(tokens[0]);
	 pos_plot.Push(tokens[1].AsInteger());

	 if(cond!="")
	 {
	    bool cond_disect=false;
	    double cond_U, cond_V, chisq;
	    if(conditionVar.Integer(tokens[0]+":"+tokens[1])==-1)
	    {
	       cond_U = SNPstat_cond.Double(SNPname);
	       cond_V = SNP_Vstat_cond.Double(SNPname);
	       chisq = cond_U*cond_U/cond_V;;
	       cond_pvalue = pchisq(chisq,1,0,0);
	       cond_effSize = cond_U/cond_V;
	       while(cond_pvalue==0.0)
	       {
		  disect=true;
		  chisq *= 0.999;
		  cond_pvalue = pchisq(chisq,1,0,0);
	       }
	    }
	    else
	    {
	       cond_effSize=0.0;
	       cond_pvalue = 1.0;
	    }
	    ifprintf(output,"%s\t%s\t%s\t%s\t%g\t%s\t%g\t%s%g\t%g\t%s%g\n",tokens[0].c_str(),tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str(),maf,direction.c_str(),effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);


	    pvalueAll_cond.Push(cond_pvalue);
	    if(maf<0.01)
	       pvalue1_cond.Push(cond_pvalue);
	    if(maf<0.05)
	       pvalue5_cond.Push(cond_pvalue);
	 }
	 else
	    ifprintf(output,"%s\t%s\t%s\t%s\t%g\t%s\t%g\t%s%g\n",tokens[0].c_str(),tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str(),maf,direction.c_str(),effSize,disect?"<":"",pvalue);

	 if(outvcf)
	    ifprintf(vcfout,"%s\t%s\t%s\t%s\t%s\t.\t.\tALT_AF=%g;\n",tokens[0].c_str(),tokens[1].c_str(),tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str(),maf);

	 //Annotate single variants to gene
	 //initialize
	 if(geneLabel.Length()==0)
	 {
	    target_chr="";
	    target_pos=0;
	    target=0;
	    target_pvalue = _NAN_;
	 }
	 if(pvalue < 0.05/SNPmaf_maf.Length() || (cond!="" && cond_pvalue < 0.05/SNPmaf_maf.Length()))
	 {
	    //if this locus is annotated already
	    if(target_chr==tokens[0] && target_pos>tokens[1].AsInteger()-1000000)
	    {
	       //if this is the higher peak then update the annotation of this locus
	       if(pvalue<target_pvalue)
	       {
		  String current_anno = group.AnnotateSingleVar(tokens[0],tokens[1].AsInteger());
		  int distance = 1000;
		  while(current_anno=="" && distance<=1000000)
		  {
		     current_anno = group.AnnotateSingleVar(tokens[0],tokens[1].AsInteger()+distance);
		     if(current_anno=="")
			current_anno = group.AnnotateSingleVar(tokens[0],tokens[1].AsInteger()-distance);
		     distance += 1000;
		  }
		  if(geneLabel[target]!="" && current_anno!="" && geneLabel[target].Find(current_anno)==-1)
		     current_anno = geneLabel[target] + "/" + current_anno;
		  if(geneLabel[target]!="" && current_anno=="")
		     current_anno = geneLabel[target];

		  geneLabel.Push(current_anno);
		  geneLabel[target]="";
		  target = geneLabel.Length()-1;
		  target_pvalue = pvalue;
		  target_pos = tokens[1].AsInteger();
		  target_chr = tokens[0];
	       }
	       //otherwise, leave the original annotation of this locus
	       else 
		  geneLabel.Push("");
	    }
	    else
	    {
	       geneLabel.Push(group.AnnotateSingleVar(tokens[0],tokens[1].AsInteger()));
	       target_chr = tokens[0];
	       target_pos = tokens[1].AsInteger();
	       target_pvalue = pvalue;
	       target = geneLabel.Length()-1;
	    }
	 }
	 else 
	    geneLabel.Push("");
      }
   }
   //Calculate genomic control
   Vector tmp;
   tmp.Copy(chisq_before_GC);
   tmp.Sort();
   double GC = tmp[0.5];
   GC /= 0.456;
   tmp.Clear();
   ifprintf(output,"#Genomic Control for pooled sample is: %g\n",GC);
   printf("  Genomic Control for all studies are:\n");
   for(int s=0;s<studies.Length();s++)
   {
      printf("  %g\t",GCbyStudy[s]);
      ifprintf(output,"#Genomic Control for study %s is: %g\n",studies[s].c_str(),GCbyStudy[s]);
   }
   //printf("\nGenomic Control for pooled sample is: %g.\n  Reminder: you can specify --correctGC option to correct for genomic control.\n",GC);

   String title,demo;
   title = "single variant analysis";
   demo = "GC=";
   demo += GC;
   writepdf.Draw(pdf,geneLabel,pvalueAll,pvalue1,pvalue5,chr_plot,pos_plot,title,demo,true);
   if(correctGC)
   {
      chisq_before_GC /= GC;
      Vector pvalue_GC_correct,pvalue1_GC,pvalue5_GC;
      pvalue_GC_correct.Dimension(chisq_before_GC.Length());
      pvalue1_GC.Dimension(pvalue1_idx.Length());
      pvalue5_GC.Dimension(pvalue5_idx.Length());
      for(int i=0;i<chisq_before_GC.Length();i++)
	 pvalue_GC_correct[i] = pchisq(chisq_before_GC[i],1,0,0);
      for(int i=0;i<pvalue1_idx.Length();i++)
	 pvalue1_GC[i] = pchisq(chisq_before_GC[pvalue1_idx[i]]/GC,1,0,0); 
      for(int i=0;i<pvalue5_idx.Length();i++)
	 pvalue5_GC[i] = pchisq(chisq_before_GC[pvalue5_idx[i]]/GC,1,0,0);
      title = "single variant analysis (GC corrected)";
      demo = "";
      writepdf.Draw(pdf,geneLabel,pvalue_GC_correct,pvalue1_GC,pvalue5_GC,chr_plot,pos_plot,title,demo,true);
   }

   if(cond!="")
   {
      title = "single variant conditional analysis";
      writepdf.Draw(pdf,geneLabel,pvalueAll_cond,pvalue1_cond,pvalue5_cond,chr_plot,pos_plot,title,demo,true);
   }
   printf("\n  done.\n\n");

   ifclose(output);
   if(outvcf) 
   {
      ifclose(vcfout);
      printf("\n  VCF file based on superset of variants from pooled studies has been saved \n    %s\n",vcf_filename.c_str());
   }
}

void Meta::Run(GroupFromAnnotation & group)
{
   if(outvcf)
      return;
   //printf("Refining SNPs in group ...\n");
   //Save maf for each annotation group to be tested for reuse
   Vector * maf;
   Vector * stats;
   Vector * cond_stats;
   Vector * singlePvalue;
   Vector * singleEffSize;
   Matrix * cov;
   Matrix * cond_cov;
   maf = new Vector [group.annoGroups.Length()];
   stats = new Vector [group.annoGroups.Length()];

   if(cond!="")
   {
      cond_stats = new Vector [group.annoGroups.Length()];
      cond_cov = new Matrix [group.annoGroups.Length()];
   }
   cov = new Matrix [group.annoGroups.Length()];
   singlePvalue = new Vector [group.annoGroups.Length()];
   singleEffSize = new Vector [group.annoGroups.Length()];

   //if(group.annoGroups.Length()<100 && !tabix) 
   //tabix = true;
   //printf("group.annoGroups.Length() is %d, tabix is %s\n",group.annoGroups.Length(),tabix?"true":"false");

   for(int g=0;g<group.annoGroups.Length();g++)
   {
      IntArray del; //del has the SNPs to be deleted from SNPlist
      maf[g].Dimension(0);
      stats[g].Dimension(0);
      int count = group.SNPlist[g].Length();
      singlePvalue[g].Dimension(0);
      singleEffSize[g].Dimension(0);

      for(int m=0;m<count;m++)
      {
	 String newSNP;
	 bool flipStatus = false;
	 double af = SNPmaf_maf[SNPmaf.Integer(group.SNPlist[g][m])];
	 double singleP = singleVarPvalue.Double(group.SNPlist[g][m]);
	 double singleEff = singleVarEff.Double(group.SNPlist[g][m]);
	 if(af==_NAN_)
	 {
	    flipStatus=true;
	    StringArray SNPflip;
	    SNPflip.AddTokens(group.SNPlist[g][m],":");
	    newSNP = SNPflip[0]+":"+SNPflip[1]+":"+SNPflip[3]+":"+SNPflip[2];
	    af = SNPmaf_maf[SNPmaf.Integer(newSNP)];
	    singleP = singleVarPvalue.Double(newSNP);
	    singleEff = singleVarEff.Double(newSNP);
	 }
	 if(af==_NAN_ || af ==0.0 || af > MAF_cutoff || singleP==_NAN_)
	 {
	    del.Push(m);
	    continue;
	 }
	 double tmp;
	 if(flipStatus)
	 {
	    group.SNPlist[g][m] = newSNP;
	    tmp = SNPstat.Double(newSNP);
	 }
	 else
	    tmp = SNPstat.Double(group.SNPlist[g][m]);
	 stats[g].Push(tmp);
	 maf[g].Push(af);

	 singlePvalue[g].Push(singleP);
	 singleEffSize[g].Push(singleEff);
      }
      //delete the SNPs that are not genotyped in any of the studies, 
      //or not polymorphic, or SNPs with maf>cutoff.
      for(int i=del.Length()-1;i>=0;i--)
      {
	 group.SNPlist[g].Delete(del[i]);
	 group.SNPNoAllele[g].Delete(del[i]);
      }
      count = group.SNPlist[g].Length();
      cov[g].Dimension(count,count,0.0);
      if(cond!="")
      {
	 cond_stats[g].Dimension(count);
	 cond_cov[g].Dimension(count,count,0.0);
      }
   }

   //loop through cov matrices of all studies and update cov
   for(int study=0;study<studies.Length();study++)
   {
      SummaryFileReader covReader;
      StringIntHash markerPosHash;
      StringArray markersInWindow,markersCov;
      String covFilename = studies[study] + ".singlevar.cov.txt.gz";
      if(FormatAdjust[study]==0)
      {
	 marker_col = 2;
	 cov_col = 3;
      }
      else
      {
	 marker_col = 4;
	 cov_col = 5;
      }
      if(tabix)
      {
	 if(!covReader.ReadTabix(covFilename))
	 {
	    covFilename = studies[study] + ".MetaCov.assoc.gz";
	    if(!covReader.ReadTabix(covFilename))
	    {
	       error("Can not open file %s. Please check file name for study %s\n",covFilename.c_str(),studies[study].c_str());
	    }
	 }
      }
      else
      {
	 printf("Reading cov matrix from study %s ...\n",studies[study].c_str());
	 String filename = studies[study] + ".singlevar.cov.txt";
	 IFILE covfile;
	 covfile  = ifopen(filename,"r");
	 if(covfile == NULL)
	 {
	    filename = studies[study] + ".singlevar.cov.txt.gz";
	    covfile  = ifopen(filename,"r");
	 }
	 if(covfile == NULL)
	 {
	    filename = studies[study] + ".MetaCov.assoc";
	    covfile  = ifopen(filename,"r");
	    marker_col = 4;
	    cov_col = 5;
	 }
	 if(covfile == NULL)
	    error("ERROR! Cannot open file: %s!\n",filename.c_str());
	 String buffer;
	 StringArray tokens;
	 int m=0;
	 while (!ifeof(covfile))
	 {
	    buffer.ReadLine(covfile);
	    if(buffer.FindChar('#')!=-1 || buffer.Find("CHROM")!=-1)
	       continue;
	    tokens.Clear();
	    tokens.AddTokens(buffer, "\t ");
	    if(tokens[0].Find("chr")!=-1)
	       tokens[0] = tokens[0].SubStr(3);
	    String SNP = tokens[0] + ":" + tokens[1];
	    //exclude variants that have different REF and ALT
	    if(SNPexclude.Integer(study+":"+SNP)!=-1)
	    {
	       continue;
	    }
	    m++;
	    markerPosHash.SetInteger(SNP,m);
	    markersInWindow.Push(tokens[marker_col]);
	    markersCov.Push(tokens[cov_col]);
	 }
	 ifclose(covfile);
	 printf("done\n");
      }
      //   printf("Updating group stats ...\n");
      //update group statistics
      for(int g=0;g<group.annoGroups.Length();g++)
      {
	 //printf("doing group %d\n",g);
	 int count = group.SNPlist[g].Length();
	 StringArray chr,pos;
	 for(int i=0;i<count;i++)
	 {
	    StringArray tmp;
	    tmp.AddTokens(group.SNPNoAllele[g][i],":");
	    chr.Push(tmp[0]);
	    pos.Push(tmp[1]);
	 } //now pos has all the positions of markers in group g.

	 Matrix cov_i,GX;
	 cov_i.Dimension(count,count,0.0);
	 if(cond!="")
	 {
	    GX.Dimension(count,XX_inv[study].cols,0.0);
	    //GX_beta.Dimension(count);
	 }
	 if(tabix)
	    for(int m=0;m<count;m++)
	    {
	       int skip = SNPexclude.Integer(study+":"+group.SNPNoAllele[g][m]);
	       if(skip!=-1)
		  //printf("SNP %s in stuy %d is skipped %d.\n",group.SNPNoAllele[g][m].c_str(),study,skip);
		  continue;
	       if(cond!="")
	       {
		  CalculateGenotypeCov(covReader,chr[m],pos[m].AsInteger(),study,GX[m]);
		  //GX_beta[m] = GX[m].InnerProduct(commonVar_betaHat[study]);
		  //cond_stats[g][m] = stats[g][m]-GX_beta[m];
		  cond_stats[g][m] = SNPstat_cond.Double(group.SNPlist[g][m]);
	       }
	       for(int s=m;s<count;s++)
	       {
		  cov_i[m][s] = GrabGenotypeCov(covReader,study,chr[m],pos[m],chr[s],pos[s],group.SNPlist[g][m],group.SNPlist[g][s]);
	       }
	    }
	 else
	    for(int m=0;m<count;m++)
	    {    
	       int loc = markerPosHash.Integer(group.SNPNoAllele[g][m]);
	       //If this SNP is not genotpyed in this study, then skip
	       if(loc==-1)
		  continue;

	       //If this SNP from this study should be excluded due to non consistent ref/alt allele
	       //then skip
	       int skip = SNPexclude.Integer(study+":"+group.SNPNoAllele[g][m]);
	       if(skip!=-1)
		  continue;

	       //if this marker is flipped then markers from the entire row 
	       //should have cov multiply by -1.0.
	       int flip = flipSNP.Integer(study+":"+group.SNPNoAllele[g][m]);
	       double multiplyFactor=1.0;
	       if(flip!=-1)
	       {
		  multiplyFactor=-1.0;
	       }
	       //read through markersInWindow and find the selected markers
	       StringArray markers,markerscov;
	       markers.AddTokens(markersInWindow[loc-1],",");
	       markerscov.AddTokens(markersCov[loc-1],",");

	       //fill in the mth row of cov
	       for(int s=m;s<count;s++)
	       {
		  int p = markers.Find(pos[s]);
		  if(p==-1)
		     continue;
		  else
		  {
		     String markerName = study+":"+chr[s]+":"+pos[s];
		     //if the marker in window is supposed to be excluded due to non-consistent ref/alt allele
		     //then skip
		     int skip = SNPexclude.Integer(markerName);
		     if(skip!=-1)
		     {
			continue;
		     }
		     int flip = flipSNP.Integer(markerName);
		     double factor=1.0;
		     if(flip!=-1)
			factor=-1.0;
		     //printf("cov=%g,N=%d\n",markerscov[p].AsDouble(),SampleSize[study]);
		     cov_i[m][s]= multiplyFactor*factor*markerscov[p].AsDouble();
		  }
	       }
	       //fill in GX
	       if(cond!="")
	       {
		  cond_stats[g][m] = SNPstat_cond.Double(group.SNPlist[g][m]);
		  String pos_str;
		  for(int s=0;s<commonVar_study[study].Length();s++)
		  {
		     if(pos[m].AsInteger()<common_pos[commonVar_study[study][s]])
			pos_str = common_pos[commonVar_study[study][s]];
		     else
		     {
			int loc = markerPosHash.Integer(common_chr[commonVar_study[study][s]]+":"+common_pos[commonVar_study[study][s]]);
			//If this SNP is not genotpyed in this study, then skip
			if(loc==-1)
			   continue;

			//read through markersInWindow and find the selected markers
			markers.Clear();
			markerscov.Clear();
			markers.AddTokens(markersInWindow[loc-1],",");
			markerscov.AddTokens(markersCov[loc-1],",");
			pos_str = pos[m];
		     }
		     int p = markers.Find(pos_str);
		     if(p==-1)
			continue;
		     GX[m][s]= markerscov[p].AsDouble();
		  }
	       }
	    }                 
	 /*
	    printf("cov_i for study %d is:\n",study);
	    for(int i=0;i<cov_i.rows;i++)
	    {
	    for(int j=0;j<cov_i.cols;j++)
	    printf("%g ",cov_i[i][j]);
	    printf("\n");
	    }
	  */
	 cov[g].Add(cov_i);
	 if(cond!="")
	 {
	    cond_cov[g].Add(cov_i);
	    Matrix GX_trans,tmp,extra_cov_i;
	    GX_trans.Transpose(GX);
	    tmp.Product(GX,XX_inv[study]);
	    extra_cov_i.Product(tmp,GX_trans);
	    extra_cov_i.Multiply(-1.0);
	    cond_cov[g].Add(extra_cov_i);
	 }
	 /*
	    for(int r=0;r<cov[g].rows;r++)
	    for(int c=r+1;c<cov[g].cols;c++)
	    {
	    cov[g][c][r] = cov[g][r][c];
	    if(cond!="")
	    cond_cov[g][c][r] = cond_cov[g][r][c];
	    }
	    printf("GX is:\n");
	    for(int i=0;i<GX.rows;i++)
	    {
	    for(int j=0;j<GX.cols;j++)
	    printf("%g\t",GX[i][j]);
	    printf("\n");
	    }
	    printf("XX_inv is:\n");
	    for(int i=0;i<XX_inv[g].rows;i++)
	    {
	    for(int j=0;j<XX_inv[g].cols;j++)
	    printf("%g\t",XX_inv[g][i][j]);
	    printf("\n");
	    }
	    printf("extra_cov_i is:\n");
	    for(int i=0;i<extra_cov_i.rows;i++)
	    {
	    for(int j=0;j<extra_cov_i.cols;j++)
	    printf("%g\t",extra_cov_i[i][j]);
	    printf("\n");
	    }
	    printf("cond_cov is:\n");
	    for(int i=0;i<extra_cov_i.rows;i++)
	    {
	    for(int j=0;j<extra_cov_i.cols;j++)
	    printf("%g\t",cond_cov[g][i][j]);
	    printf("\n");
	    }
	    printf("cov_i is:\n");
	    for(int i=0;i<extra_cov_i.rows;i++)
	    {
	    for(int j=0;j<cov_i.cols;j++)
	    printf("%g\t",cov_i[i][j]);
	    printf("\n");
	    }
	    printf("cov is:\n");
	    for(int i=0;i<cov[g].rows;i++)
	    {
	    for(int j=0;j<cov[g].cols;j++)
	    printf("%g\t",cov[g][i][j]);
	    printf("\n");
	    }
	    printf("stats and cond_stats are:\n");
	    for(int i=0;i<cond_stats[g].Length();i++)
	    printf("%g %g\t",stats[g][i],cond_stats[g][i]);
	    printf("\n");
	  */
      }
   }

   String method = "";
   if(Burden)
   {
      method = "burden";
      BurdenAssoc(method,group,maf,stats,cond_stats,cov,cond_cov,singleEffSize,singlePvalue);
   }
   if(MB)
   {
      method = "MB";
      BurdenAssoc(method,group,maf,stats,cond_stats,cov,cond_cov,singleEffSize,singlePvalue);
   }

   if(VTa)
   {
      printf("Performing Variable Threshold tests ...\n");
      //calculate final results here
      Vector pvalue_VT,pos_plot;
      StringArray chr_plot,geneLabels;

      IFILE output;
      String filename;
      if(prefix =="")
	 filename = "meta.VT.results";
      else if(prefix.Last()=='.' || prefix.Last()=='/')
	 filename = prefix +  "meta.VT.results";
      else
	 filename = prefix + ".meta.VT.results";
      output=ifopen(filename,"w",InputFile::UNCOMPRESSED);

      String method = "VT_";
      method += MAF_cutoff;
      IFILE reportOutput;
      if(report)
      {
	 String reportFile;
	 if(prefix =="")
	    reportFile = "meta.tophits.VT.tbl";
	 else if(prefix.Last()=='.' || prefix.Last()=='/')
	    reportFile = prefix +  "meta.tophits.VT.tbl";
	 else
	    reportFile = prefix + ".meta.tophits.VT.tbl";
	 reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
	 ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
      }

      ifprintf(output,"##Method=VT\n");
      ifprintf(output,"##STUDY_NUM=%d\n",studies.Length());
      ifprintf(output,"##TotalSampleSize=%d\n",total_N);
      if(fullResult)
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\n");
      else
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\n");

      for(int g=0;g<group.annoGroups.Length();g++)
      {
	 if(g>1 && g%1000==1)
	    printf("Finished analyzing %d genes.\n",g-1);

	 if(maf[g].Length()==0)
	 {
	    if(fullResult)
	       ifprintf(output,"%s\t0\t-\t-\t-\t-\t0\t0\t0\t-\t-\t-\n",group.annoGroups[g].c_str());
	    else
	       ifprintf(output,"%s\t0\t-\t0\t0\t0\t-\t-\t-\n",group.annoGroups[g].c_str());
	    continue;
	 }

	 if(maf[g].Length()==1)
	 {
	    if(fullResult)
	    {
	       ifprintf(output,"%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),group.SNPlist[g][0].c_str(),maf[g][0],singleEffSize[g][0],singlePvalue[g][0],maf[g][0],maf[g][0],maf[g][0],singleEffSize[g][0],maf[g][0],singlePvalue[g][0]);
	    }
	    else
	    {
	       ifprintf(output,"%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),group.SNPlist[g][0].c_str(),maf[g][0],maf[g][0],maf[g][0],singleEffSize[g][0],maf[g][0],singlePvalue[g][0]);
	    }
	    pvalue_VT.Push(singlePvalue[g][0]);

	    StringArray tmp_SNPname;
	    tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
	    chr_plot.Push(tmp_SNPname[0]);
	    pos_plot.Push(tmp_SNPname[1].AsDouble());
	    geneLabels.Push(group.annoGroups[g]);
	    continue;
	 }

	 //STEP1: sort maf[g] to find # of cutoffs
	 Vector cp_maf;
	 cp_maf.Copy(maf[g]);
	 cp_maf.Sort();
	 Vector maf_cutoff;
	 maf_cutoff.Push(cp_maf[0]);

	 for(int i=1;i<cp_maf.Length();i++)
	 {
	    if(cp_maf[i]>maf_cutoff[maf_cutoff.Length()-1])
	    {
	       maf_cutoff.Push(cp_maf[i]);	
	    }
	 } //now unique maf cutoffs are saved in maf_cutoff.

	 //STEP2: do burden test for each grouping method and find t_max
	 double pvalue,chosen_cutoff,chosen_effSize;
	 double numerator,denominator,t_max=_NAN_;
	 Vector weight,tmp,chosen_weight,score;
	 Matrix cov_weight;
	 weight.Dimension(maf[g].Length());
	 tmp.Dimension(group.SNPlist[g].Length());
	 cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());
	 score.Dimension(maf_cutoff.Length());

	 for(int i=0;i<maf_cutoff.Length();i++)
	 {
	    for(int w=0;w<weight.Length();w++)
	    {
	       if(maf[g][w]<=maf_cutoff[i])
		  weight[w]=1.0;
	       else
		  weight[w]=0.0;
	       cov_weight[i][w]=weight[w];
	    }
	    numerator = weight.InnerProduct(stats[g]);

	    for(int d=0;d<tmp.Length();d++)
	    {
	       tmp[d] = weight.InnerProduct(cov[g][d]);
	    }
	    denominator = tmp.InnerProduct(weight);

	    if(denominator != 0)
	    {
	       double t_stat = fabs(numerator/sqrt(denominator));
	       score[i]=t_stat;
	       if(t_max==_NAN_)
	       {
		  t_max = t_stat;
		  chosen_cutoff = maf_cutoff[i];
		  chosen_weight.Copy(weight);
		  chosen_effSize = numerator/denominator;
	       }
	       else
	       {
		  if(t_stat>t_max) 
		  {
		     t_max = t_stat;
		     chosen_cutoff = maf_cutoff[i];
		     chosen_weight.Copy(weight);
		     chosen_effSize = numerator/denominator;
		  }
	       }
	    }
	    else 
	       score[i]=0.0;
	 }

	 if(score.Max()==0.0)
	 {
	    if(fullResult)
	       ifprintf(output,"%s\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n",group.annoGroups[g].c_str());
	    else
	       ifprintf(output,"%s\t-\t-\t-\t-\t-\t-\t-\t-\n",group.annoGroups[g].c_str());
	    continue;
	 }
	 Vector tmp_maf,tmp_eff,tmp_pvalue;
	 for(int i=0;i<maf[g].Length();i++)
	 {
	    if(chosen_weight[i]==1.0)
	    {
	       tmp_maf.Push(maf[g][i]);
	    }
	 }

	 for(int i=0;i<maf[g].Length();i++)
	 {
	    if(chosen_weight[i]==1.0)
	    {
	       tmp_eff.Push(singleEffSize[g][i]);
	       tmp_pvalue.Push(singlePvalue[g][i]);
	    }
	 }

	 double average_af = tmp_maf.Average();
	 double min_af = tmp_maf.Min();
	 double max_af = tmp_maf.Max();

	 String var;
	 for(int i=0;i<maf[g].Length()-1;i++)
	 {
	    if(chosen_weight[i]==1.0)
	       var += group.SNPlist[g][i] + ";";
	 }
	 if(chosen_weight[maf[g].Length()-1]==1.0)
	    var += group.SNPlist[g][maf[g].Length()-1];

	 //STEP3: calculate covariance matrix for (U_1 ... U_#cutoff) 
	 Matrix cov_U,cov_U_tmp;
	 cov_U_tmp.Product(cov_weight,cov[g]);
	 Matrix cov_weight_trans;
	 cov_weight_trans.Transpose(cov_weight);
	 cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
	 //Calculate covariance matrix for (T_1 ... T_#cutoff)
	 Matrix cov_T;
	 cov_T.Dimension(cov_U.rows,cov_U.cols);
	 cov2cor(cov_U,cov_T);

	 //STEP4: calculate VT pvalue and report.
	 int cutoff = maf_cutoff.Length();
	 double * lower = new double [cutoff];
	 double * upper = new double [cutoff];
	 double * mean = new double [cutoff];

	 for(int i=0;i<cutoff;i++)
	 {
	    mean[i] = 0.0;
	    lower[i] = -t_max;
	    upper[i] = t_max;
	 }

	 //Use pmvnorm to calculate the asymptotic p-value
	 Vector result;
	 pmvnorm(lower,upper,mean,cov_T,false,result);
	 //if(result[0]==-1.0)
	 if(result[0]==-1.0)
	 {
	    if(fullResult)
	    {
	       ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str());

	       for(int i=0;i<tmp_maf.Length()-1;i++)
		  ifprintf(output,"%g,",tmp_maf[i]);
	       ifprintf(output,"%g\t",tmp_maf[tmp_maf.Length()-1]);

	       for(int i=0;i<tmp_eff.Length()-1;i++)
		  ifprintf(output,"%g,",tmp_eff[i]);
	       ifprintf(output,"%g\t",tmp_eff[tmp_eff.Length()-1]);

	       for(int i=0;i<tmp_pvalue.Length()-1;i++)
		  ifprintf(output,"%g,",tmp_pvalue[i]);
	       ifprintf(output,"%g\t",tmp_pvalue[tmp_pvalue.Length()-1]);

	       ifprintf(output,"%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\n",average_af,min_af,max_af,chosen_effSize,chosen_cutoff);
	    }
	    else
	       ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\n",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str(),average_af,min_af,max_af,chosen_effSize,chosen_cutoff);
	 }
	 else
	 {
	    if(1.0-result[0]==0.0)
	    {
	       //	    printf("gene %s has result %g\n",group.annoGroups[g].c_str(),1.0-result[0]);
	       printf("Using Shuang's algorithm to calculate MVN pvalue for gene %s ... ",group.annoGroups[g].c_str());
	       if(maf_cutoff.Length()>20)
	       {
		  while(maf_cutoff.Length()>20)
		  {
		     maf_cutoff.Delete(0);
		  }

		  double numerator,denominator,t_max=_NAN_;
		  Vector weight,tmp,chosen_weight,score;
		  Matrix cov_weight;
		  weight.Dimension(maf[g].Length());
		  tmp.Dimension(group.SNPlist[g].Length());
		  cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());

		  for(int i=0;i<maf_cutoff.Length();i++)
		  {
		     for(int w=0;w<weight.Length();w++)
		     {
			if(maf[g][w]<=maf_cutoff[i])
			   weight[w]=1.0;
			else
			   weight[w]=0.0;
			cov_weight[i][w]=weight[w];
		     }
		     numerator = weight.InnerProduct(stats[g]);

		     for(int d=0;d<tmp.Length();d++)
		     {
			tmp[d] = weight.InnerProduct(cov[g][d]);
		     }
		     denominator = tmp.InnerProduct(weight);
		     if(denominator != 0)
		     {
			double t_stat = fabs(numerator/sqrt(denominator));
			score.Push(t_stat);
			if(t_max==_NAN_)
			{
			   t_max = t_stat;
			   chosen_cutoff = maf_cutoff[i];
			   chosen_weight.Copy(weight);
			   chosen_effSize = numerator/denominator;
			}
			else
			{
			   if(t_stat>t_max)
			   {
			      t_max = t_stat;
			      chosen_cutoff = maf_cutoff[i];
			      chosen_weight.Copy(weight);
			      chosen_effSize = numerator/denominator;
			   }
			}
		     }
		     else
			score.Push(0.0);
		  }
		  if(score.Max()==0.0)
		  {
		     if(fullResult)
			ifprintf(output,"%s\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n",group.annoGroups[g].c_str());
		     else
			ifprintf(output,"%s\t-\t-\t-\t-\t-\t-\t-\t-\n",group.annoGroups[g].c_str());
		     continue;
		     printf("completed!\n");
		  }

		  Vector tmp_maf,tmp_eff,tmp_pvalue;
		  for(int i=0;i<maf[g].Length();i++)
		  {
		     if(chosen_weight[i]==1.0)
		     {
			tmp_maf.Push(maf[g][i]);
		     }
		  }

		  for(int i=0;i<maf[g].Length();i++)
		  {
		     if(chosen_weight[i]==1.0)
		     {
			tmp_eff.Push(singleEffSize[g][i]);
			tmp_pvalue.Push(singlePvalue[g][i]);
		     }
		  }
		  average_af = tmp_maf.Average();
		  min_af = tmp_maf.Min();
		  max_af = tmp_maf.Max();

		  String var;
		  for(int i=0;i<maf[g].Length()-1;i++)
		  {
		     if(chosen_weight[i]==1.0)
			var += group.SNPlist[g][i] + ";";
		  }
		  if(chosen_weight[maf[g].Length()-1]==1.0)
		     var += group.SNPlist[g][maf[g].Length()-1];

		  //STEP3: calculate covariance matrix for (U_1 ... U_#cutoff) 
		  Matrix cov_U,cov_U_tmp;
		  cov_U_tmp.Product(cov_weight,cov[g]);
		  Matrix cov_weight_trans;
		  cov_weight_trans.Transpose(cov_weight);
		  cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
		  //Calculate covariance matrix for (T_1 ... T_#cutoff)
		  Matrix cov_T;
		  cov_T.Dimension(cov_U.rows,cov_U.cols);
		  cov2cor(cov_U,cov_T);

		  pvalue = CalculateMVTPvalue(score,cov_T,t_max);
		  printf("completed!\n");
	       }
	       else
	       {
		  pvalue = CalculateMVTPvalue(score,cov_T,t_max);
		  printf("completed!\n");
	       }
	    }
	    else
	    {
	       pvalue = 1.0 - result[0];
	    }

	    if(pvalue <report_pvalue_cutoff && report)
	    {
	       StringArray variants;
	       variants.AddTokens(var,";");
	       for(int v=0;v<tmp_maf.Length();v++)
	       {
		  ifprintf(reportOutput,"%s\t%s\t%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),pvalue,MAF_cutoff,chosen_cutoff,variants[v].c_str(),tmp_maf[v],tmp_eff[v],tmp_pvalue[v]);
	       }
	    }

	    if(fullResult)
	    {
	       ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str());

	       for(int i=0;i<tmp_maf.Length()-1;i++)
		  ifprintf(output,"%g,",tmp_maf[i]);
	       ifprintf(output,"%g\t",tmp_maf[tmp_maf.Length()-1]);

	       for(int i=0;i<tmp_eff.Length()-1;i++)
		  ifprintf(output,"%g,",tmp_eff[i]);
	       ifprintf(output,"%g\t",tmp_eff[tmp_eff.Length()-1]);

	       for(int i=0;i<tmp_pvalue.Length()-1;i++)
		  ifprintf(output,"%g,",tmp_pvalue[i]);
	       ifprintf(output,"%g\t",tmp_pvalue[tmp_pvalue.Length()-1]);

	       ifprintf(output,"%g\t%g\t%g\t%g\t%g\t%g\n",average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
	    }
	    else
	       ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str(),average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
	    if(pvalue>1.0) 
	       pvalue = 1.0;
	 }

	 pvalue_VT.Push(pvalue);

	 StringArray tmp_SNPname;
	 tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
	 chr_plot.Push(tmp_SNPname[0]);
	 pos_plot.Push(tmp_SNPname[1].AsDouble());
	 geneLabels.Push(group.annoGroups[g]);

	 if(lower) delete [] lower;
	 if(upper) delete [] upper;
	 if(mean) delete [] mean;
      }

      String name = "VT (maf<";
      name +=  MAF_cutoff;
      name +=  ")";
      String extraname = "";
      String demo="";
      writepdf.Draw(pdf,geneLabels,pvalue_VT,chr_plot,pos_plot,name,extraname,demo,true);

      ifclose(output);
      if(report)
	 ifclose(reportOutput);
      printf("Done.\n\n");
   }

   if(SKAT)
   {
      printf("Performing SKAT ...\n");
      //calculate Q statistics here

      double *rho = new double [8];
      if(SKATO)
      {
	 rho[0]=0.0; rho[1]=0.01; rho[2]=0.04;rho[3]=0.9;
	 rho[4]=0.16; rho[5]=0.25; rho[6]=0.5; rho[7]=1.0;
      }
      Vector pvalue_SKAT,pos_plot,cond_pvalue_SKAT;
      StringArray chr_plot,geneLabels;
      IFILE output;
      String filename;
      if(prefix =="")
	 filename = "meta.SKAT.results";
      else if(prefix.Last()=='.' || prefix.Last()=='/')
	 filename = prefix +  "meta.SKAT.results";
      else
	 filename = prefix + ".meta.SKAT.results";
      output=ifopen(filename,"w",InputFile::UNCOMPRESSED);

      String method = "SKAT_";
      method += MAF_cutoff;
      IFILE reportOutput;
      if(report)
      {
	 String reportFile;
	 if(prefix =="")
	    reportFile = "meta.tophits.SKAT.tbl";
	 else if(prefix.Last()=='.' || prefix.Last()=='/')
	    reportFile = prefix +  "meta.tophits.SKAT.tbl";
	 else
	    reportFile = prefix + ".meta.tophits.SKAT.tbl";
	 reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
	 ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE_DAVIES\tGENE_PVALUE_LIU\tMAF_CUTOFF\tACTUAL_CUTOFF\tVAR\tMAF\tEFFSIZE\tPVALUE\n");
      }

      ifprintf(output,"##Method=SKAT\n");
      ifprintf(output,"##STUDY_NUM=%d\n",studies.Length());
      ifprintf(output,"##TotalSampleSize=%d\n",total_N);
      if(fullResult)
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");
      else
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tSTATISTICS\tPVALUE_DAVIES\tPVALUE_LIU\n");

      double Qstat,pvalue,pvalue_liu;
      for(int g=0;g<group.annoGroups.Length();g++)
      {
	 //      printf("Now working on group %s\n",group.annoGroups[g].c_str());
	 if(g%1000==1 && g>1000)
	 {
	    printf("Finished analyzing %d genes.\n",((int)g/1000)*1000);
	 }
	 Vector weight;
	 int n = maf[g].Length();
	 if(maf[g].Length()==0)
	 {
	    if(fullResult)
	       ifprintf(output,"%s\t0\t-\t-\t-\t-\t0\t0\t0\t-\t-\t-\n",group.annoGroups[g].c_str());
	    else
	       ifprintf(output,"%s\t0\t-\t0\t0\t0\t-\t-\t-\n",group.annoGroups[g].c_str());
	    continue;
	 }

	 double average_af = maf[g].Average();
	 double min_af = maf[g].Min();
	 double max_af = maf[g].Max();

	 String var;
	 for(int i=0;i<maf[g].Length()-1;i++)
	 {
	    var += group.SNPlist[g][i] + ";";
	 }
	 var += group.SNPlist[g][maf[g].Length()-1];

	 //refine variants  
	 Vector cp_maf;
	 cp_maf.Copy(maf[g]);

	 if(n>1)
	 {
	    //calculate correlation of cov[g]_i and cov[g]_j
	    IntArray del;
	    double corr;
	    for(int i=0;i<n;i++)
	    {
	       for(int j=i+1;j<n;j++)
	       {
		  corr = CalculateCorrCoef(cov[g][i],cov[g][j]);
		  if(corr ==1.0 || corr == -1.0)
		  {
		     del.Push(i);
		     break;
		  }
	       }
	    }
	    //del.RemoveDuplicates();

	    for(int i=del.Length()-1;i>=0;i--)
	    {
	       cov[g].DeleteRow(del[i]);
	       cov[g].DeleteColumn(del[i]);
	       stats[g].Delete(del[i]);
	       cp_maf.Delete(del[i]);
	       if(cond!="")
	       {
		  cond_stats[g].Delete(del[i]);
		  cond_cov[g].DeleteRow(del[i]);
		  cond_cov[g].DeleteColumn(del[i]);
	       }
	       n--;
	    }
	 }
	 weight.Dimension(n);
	 //get weight based on maf
	 double alpha = 1.0;
	 double beta=25.0;
	 double beta_density;
	 for(int w=0;w<n;w++)
	 {
	    beta_density = GetBetaDensity(alpha,beta,cp_maf[w]);
	    weight[w] = beta_density * beta_density;
	 }
	 Vector tmp,cond_tmp;
	 tmp.Dimension(n);
	 if(cond!="")
	    cond_tmp.Dimension(n);
	 for(int i=0;i<n;i++)
	 {
	    tmp[i] = weight[i]*stats[g][i];
	    if(cond!="")
	       cond_tmp[i] = weight[i]*cond_stats[g][i];
	 }
	 Qstat = tmp.InnerProduct(stats[g]);
	 /*
	    printf("Qstat is: %g\n U and w are:\n",Qstat);
	    for(int i=0;i<stats[g].dim;i++)
	    {
	    printf("%g\t%g\n",stats[g][i],weight[i]);
	    }
	    printf("\n");
	  */
	 double cond_Qstat;
	 if(cond!="")
	    cond_Qstat = cond_tmp.InnerProduct(cond_stats[g]);
	 double * lambda = new double [n];
	 CalculateLambda(cov[g],weight,lambda,false,0.0);
	 /*
	    printf("lambdas are :\n");
	    for(int i=0;i<n;i++)
	    printf("%g\t",lambda[i]);
	    printf("\n");
	  */
	 double Qstat_dav = Qstat;
	 double Qstat_liu = Qstat;
	 double cond_pvalue,cond_pvalue_liu;

	 pvalue = MixChidist(lambda, n, Qstat,"Davies");

	 bool disect_davies=false,disect_liu = false,cond_disect_davies=false,cond_disect_liu=false;
	 int disect_itr=0;
	 while( (pvalue<=0.0 ||pvalue==2.0 || isnan(pvalue)) && disect_itr<10000)
	 {
	    disect_davies=true;
	    Qstat_dav*=0.9999;
	    pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
	    disect_itr++;
	 }
	 while((pvalue<=0.0 ||pvalue==2.0 || isnan(pvalue)))
	 {
	    Qstat_dav*=0.99;
	    pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
	 }
	 pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	 disect_itr=0;
	 while( (pvalue_liu<=0.0 ||pvalue_liu==2.0 || isnan(pvalue_liu)) && disect_itr<10000)
	 {
	    disect_liu=true;
	    Qstat_liu*=0.9999;
	    pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	    disect_itr++;
	 }
	 while((pvalue_liu<=0.0 ||pvalue_liu==2.0 || isnan(pvalue_liu)))
	 {
	    Qstat_liu*=0.99;
	    pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	 }

	 if(cond!="")
	 {
	    double * lambda = new double [n];
	    CalculateLambda(cond_cov[g],weight,lambda,false,0.0);
	    Qstat_dav = cond_Qstat;
	    Qstat_liu = cond_Qstat;

	    cond_pvalue = MixChidist(lambda, n, cond_Qstat,"Davies");

	    int disect_itr=0;
	    while( (cond_pvalue<=0.0 ||cond_pvalue==2.0 || isnan(cond_pvalue)) && disect_itr<10000)
	    {
	       cond_disect_davies=true;
	       Qstat_dav*=0.9999;
	       pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
	       disect_itr++;
	    }
	    while((cond_pvalue<=0.0 ||cond_pvalue==2.0 || isnan(cond_pvalue)))
	    {
	       Qstat_dav*=0.99;
	       cond_pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
	    }
	    cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	    disect_itr=0;
	    while( (cond_pvalue_liu<=0.0 ||cond_pvalue_liu==2.0 || isnan(cond_pvalue_liu)) && disect_itr<10000)
	    {
	       cond_disect_liu=true;
	       Qstat_liu*=0.9999;
	       cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	       disect_itr++;
	    }
	    while((cond_pvalue_liu<=0.0 ||cond_pvalue_liu==2.0 || isnan(cond_pvalue_liu)))
	    {
	       Qstat_liu*=0.99;
	       cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	    }
	    if(lambda) delete [] lambda;
	 }
	 if(isnan(pvalue_liu) || isnan(pvalue))
	 {
	    if(fullResult)
	    {
	       ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());

	       for(int i=0;i<maf[g].Length()-1;i++)
		  ifprintf(output,"%g,",maf[g][i]);
	       ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);

	       for(int i=0;i<singleEffSize[g].Length()-1;i++)
		  ifprintf(output,"%g,",singleEffSize[g][i]);
	       ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);

	       for(int i=0;i<singlePvalue[g].Length()-1;i++)
		  ifprintf(output,"%g,",singlePvalue[g][i]);
	       ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
	       if(cond=="")
		  ifprintf(output,"%g\t%g\t%g\t%g\t-\t-\n",average_af,min_af,max_af,Qstat);
	       else
		  ifprintf(output,"%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n",average_af,min_af,max_af,Qstat,cond_Qstat);
	    }
	    else
	       if(cond=="")
		  ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat);
	       else
		  ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t-\t-\t%g\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,cond_Qstat);
	    continue;
	 }
	 //tabulate top results
	 if(cond!="")
	    if((cond_pvalue <report_pvalue_cutoff && cond_pvalue_liu<report_pvalue_cutoff)  && report)
	    {
	       StringArray variants;
	       variants.AddTokens(var,";");
	       for(int v=0;v<maf[g].Length();v++)
	       {
		  ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
	       }
	    }
	    else
	       if((pvalue <report_pvalue_cutoff && pvalue_liu<report_pvalue_cutoff)  && report)
	       {
		  StringArray variants;
		  variants.AddTokens(var,";");
		  for(int v=0;v<maf[g].Length();v++)
		  {
		     ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
		  }
	       }

	 if(fullResult)
	 {
	    ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());

	    for(int i=0;i<maf[g].Length()-1;i++)
	       ifprintf(output,"%g,",maf[g][i]);
	    ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);

	    for(int i=0;i<singleEffSize[g].Length()-1;i++)
	       ifprintf(output,"%g,",singleEffSize[g][i]);
	    ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);

	    for(int i=0;i<singlePvalue[g].Length()-1;i++)
	       ifprintf(output,"%g,",singlePvalue[g][i]);
	    ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
	    if(cond=="")
	       ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%s%g\n",average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu);
	    else
	       ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n",average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,cond_Qstat,cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu);
	 }
	 else
	 {
	    if(cond=="")
	       ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu);
	    else
	       ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%s%g\t%g\t%s%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,Qstat,disect_davies?"<":"",pvalue,disect_liu?"<":"",pvalue_liu,cond_Qstat,cond_disect_davies?"<":"",cond_pvalue,cond_disect_liu?"<":"",cond_pvalue_liu);
	 }

	 if(lambda) delete [] lambda;

	 pvalue_SKAT.Push(pvalue);
	 if(cond!="")
	    cond_pvalue_SKAT.Push(cond_pvalue);
	 StringArray tmp_SNPname;
	 tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
	 chr_plot.Push(tmp_SNPname[0]);
	 pos_plot.Push(tmp_SNPname[1].AsDouble());
	 geneLabels.Push(group.annoGroups[g]);
      }

      String name = "SKAT (maf<";
      name +=  MAF_cutoff;
      name +=  ")";
      String extraname = "";
      String demo = "";
      writepdf.Draw(pdf,geneLabels,pvalue_SKAT,chr_plot,pos_plot,name,extraname,demo,true);
      if(cond!="")
      {
	 name += " conditional analysis";
	 writepdf.Draw(pdf,geneLabels,cond_pvalue_SKAT,chr_plot,pos_plot,name,extraname,demo,true);
      }
      ifclose(output);
      if(report)
	 ifclose(reportOutput);
      printf("Done.\n\n");
      if(rho) delete [] rho;
   }

   pdf.CloseFile();

   //housekeeping
   if(maf) delete [] maf;
   if(stats) delete [] stats;
   if(cov) delete [] cov;
   if(singlePvalue) delete [] singlePvalue;
   if(singleEffSize) delete [] singleEffSize;
   printf("\nQQ plots and manhattan polts have been saved in %s.\n",pdf_filename.c_str());
}

double Meta::CalculateSKATOpvalue(double * rho,Matrix & cov,Vector & weight)
{
   int n=cov.rows;
   //n.q<-dim(Q.all)[1]
   //p.m<-dim(Phi)[2]
   //calcualte Lambda for each rho
   double ** lambda = new double * [5];
   dobule * lambda = new double [n];
   //construct Phi from cov and weight
   Matrix Phi;
Phi.Dimension(cov.rows,cov.cols);

   Vector sqrt_w(weight.Length());
   for(int i=0;i<cov.rows;i++)
      sqrt_w[i] = sqrt(weight[i]);
   for(int i=0;i<cov.rows;i++)
      for(int j=0;j<cov.cols;j++)
	 Phi[i][j] = sqrt_w[i]*cov[i][j]*sqrt_w[j];

   Cholesky chol;
   for(int i=0;i<7;i++)
   {
      //construct R(rho)
      Matrix R,phi_rho,tmp;
      R.Dimension(n,n,rho[i]);
      for(int i=0;i<n;i++)
	 R[i][i] =1.0;
      //cholesky decomp of R
      FastDecompose(R);
      //calculate phi(rho)
      tmp.Product(chol.L,Phi);
      Matrix transL;
      transL.Transpose(chol.L);
      phi_rho.Product(tmp,transL);
      //Calculate lambda for phi_rho
      CalculateSKATOLambda(phi_rho,lambda[i]);
   }
   double MuQ,VarQ,KerQ,VarRemain,Df,tau,tau1,tau2,z_mean_2;
   double lambdaOpt[Phi.rows];
   int m;
   double * Q_all;
   for(int i=0;i<7;i++)
{
      Q_all[i] = (1.0-rho[i])*Q_SKAT + rho[i] * Q_burden;
Q_all[i] /= 2.0;
 }
  //Calculate parameters for SKATO
   CalculateSKATOParam(Phi,rho,MuQ,VarQ,KerQ,lambdaOpt,VarRemain,Df,tau,z_mean_2,m, tau1,tau2);
   SKAT_Optimal_Each_Q(m,Q_all,rho,lambda);
   //finalize SKATO pvalue
   //integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=40, subdivisions=1000, pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-25), silent = TRUE;

   double a, b, epsabs, epsrel, result, abserr, *work;
   int neval, ier, last, *iwork;

   // a and b are the integration limits
   a = 0.0;
   b = 40;

   // epsabs and epsrel determine the accuracy requirement
   epsabs = 10e-15;
   epsrel = 10e-15;

   // allocate memory
   iwork = (int *) calloc(LIMIT, sizeof(int));
   work = (double *) calloc(LENW, sizeof(double));

   double result = cquad_dqags(FuncToIntegrate, a, b, epsabs, epsrel, &abserr, &neval, &ier, &last, iwork, work);
   return 1.0-result;
}

void Meta::CalculateSKATOLambda(Matrix & phi_rho,double * lambda)
{
   SVD svd;
   svd.Decompose(phi_rho);
   for(int i=0;i<phi_rho.rows;i++)
      lambda[i] = svd.w[i];
}

void Meta::CalculateSKATOParam(Matrix & Phi,double * rho,double & MuQ,double & VarQ, double & KerQ,double * lambda,double & VarRemain,double & Df,double & tau,double & z_mean_2,double &m,double & tau1,double & tau2)
{
   m = Phi.cols;
   int n = 7;
   // ZMZ
   double denom=0.0;
   for(int i=0;i<m;i++)
      for(int j=0;j<m;j++)
	 denom += Phi[i][j];
   double * Z1_1;
   for(int i=0;i<m;i++)
      Z1_1[i] = 0.0;
   for(int i=0;i<m;i++)
   {
      for(int j=0;j<m;j++)
	 Z1_1[i] += Phi[i][i];
   }
   Matrix ZMZ;
   ZMZ.Dimension(m,m);
   for(int i=0;i<m;i++)
      for(int j=0;j<m;j++)
	 ZMZ[i][j] = Z1_1[i]*Z1_1[j]/denom;

   // W3.2 Term : mixture chisq
   Matrix W3_2_t,tmp;
   W3_2_t = Phi - ZMZ;
   tmp.Product(ZMZ,W3_2-t);
   CalculateSKATOLambda(W3_2_t,lambda);

   //W3.3 Term : variance of remaining ...
   VarRemain=0.0;
   for(int i=0;i<m;i++)
      for(int j=0;j<m;j++)
	 VarRemain += tmp[i][j];
   VarRemain *= 4.0;

   //tau term
   z_mean_2 = denom/(m*m);
   tmp.Product(Phi,Phi);
   double sum=0.0;
   for(int i=0;i<m;i++) 
      for(int j=0;j<m;j++) 
	 sum += tmp[i][j]; 
   tau1 =  sum / (m*m) / z_mean_2;

   // Mixture Parameters
   MuQ=0.0;
   double sum_sqr=0.0, sum_quad=0.0;
   for(int i=0;i<m;i++)
   {
      MuQ += lambda[i]; 
      double sqr = lambda[i]*lambda[i];
      sum_sqr += sqr;
      sqr *= sqr;
      sum_quad += sqr;
   }
   VarQ = sum_sqr*2.0 + VarRemain;
   KerQ = sum_quad*12.0/(sum_sqr*sum_sqr);
   Df = 12/KerQ;

   for(int i=0;i<n;i++)
      tau[i]=m*m*rho[i] * z_mean_2 + tau1 * (1.0-rho[i]);
   tau2 = m*z_mean_2;
}
void Meta::Get_Liu_Params_Mod(double *c1,double & muQ,double & varQ,double & df)
{
   muQ = c1[0];
   varQ = 2.0 *c1[1];
   double s1 = c1[2] / pow(c1[1],(3.0/2.0));
   double s2 = c1[3] / c1[1]*c1[1] ;

   double beta1 = sqrt(8.0)*s1;
   double beta2 = 12.0*s2;
   double type1 = 0.0;
   double a,d;

   if(pow(s1,2) > s2)
   {
      a = 1.0/(s1 - sqrt(s1*s1 - s2));
      d = s1 *a*a*a - a*a;
      df = a*a - 2.0*d;
   } else {
      type1 = 1;
      df = 1/s2;
      a = sqrt(l);
      d = 0;
   }
   //double muX  = df+d;
   //double sigmaX = sqrt(2.0) *a;
}

void Meta::CalculateOptimalEachQ(int m,double * Q_all, double ** lambda,double &pmin,double * qmin)
{
   int n.r = 7,n.q=7;

   double pval[7];
   double qmin[7];
   double param[7][3];

   for(int i=0;i<7;i++)
   {
      double c1[4];
      c1[0]=0.0;c1[1]=0.0;c1[2]=0.0;c1[3]=0.0;
      for(int j=0;j<m;j++)
      {
	 c1+= lambda[i][j];
	 double sqr =  lambda[i][j]*lambda[i][j];
	 c1[1] += sqr;
	 c1[2] += sqr*lambda[i][j];
	 c1[3] += sqr*sqr;
      }
      Get_Liu_Params_Mod(c1,param[i][0],param[i][1],param[i][2]);
      pval[i] =  MixChidist(lambda[i], m, Q_all[i],"Davies");
   }
   pmin = pval[0];
   for(int i=0;i<7;i++)
   {
      if(pval[i]<pmin)
	 pmin=pval[i];
   }
   for(int i=0;i<7;i++)
   {

      double q.org = qchisq(1.0-pmin,param[i][2]);
      qmin[i] = (q.org - param[i][2])/sqrt(2.0*param[i][2]) *sqrt(param[i][1]) + param[i][0];
   }
   if(param) delete []param;
}

double Meta::FuncToIntegrate(double *x,double *rho,double *qmin,double *tau,double &muQ,double &VarQ,double &VarRemain)
{
   double result,q.org,q.q;
   double min = 6.66e66;
   double temp;
   for(int i=0;i<7;i++)
   {
      temp = (qmin[i] - *x * tau[i])/(1.0-rho[i]);
      if(temp<min)
	 min=temp;
   }
   double sd = sqrt((VarQ-VarRemain)/VarQ);
   min -= muQ;
   min *=sd;
   min += muQ;
   //calculate pvalue
   result =  MixChidist(lambda, n, min,"Davies");
   if(result>1.0)
      result =1.0;
   result = (1.0-result)*dchisq(*x,1);
   return result;
}

void Meta::CalculateLambda(Matrix & cov,Vector & weight, double * lambda,bool optimal,double rho)
{
   SVD svd;
   int n = cov.rows;
   Matrix L,final;
   L.Dimension(n,n);
   final.Dimension(n,n);
   //calculat sqrt(V)
   svd.Decompose(cov);
   Matrix tmp,tmp2;
   tmp.Dimension(n,n);
   tmp2.Dimension(n,n);
   for(int i=0;i<n;i++)
   {
      for(int j=0;j<n;j++)
	 tmp[i][j] = svd.u[i][j]*sqrt(svd.w[j]);
   }
   Matrix v;
   v.Transpose(svd.u);
   tmp2.Product(tmp,v);

   if(!optimal)
   {
      for(int i=0;i<n;i++)
      {
	 for(int j=0;j<n;j++)
	    tmp[i][j] = tmp2[i][j] * weight[j];
      }
      final.Product(tmp,tmp2);
      svd.Decompose(final);
      for(int i=0;i<n;i++)
	 lambda[i] = fabs(svd.w[i]);
   } 
   else
   {
      Matrix R;
      R.Dimension(n,n,rho);
      for(int i=0;i<n;i++)
	 R[i][i] = 1.0;

      for(int i=0;i<n;i++)
      {
	 for(int j=0;j<n;j++)
	    tmp[i][j] = tmp2[i][j] * sqrt(weight[j]);
      }
      Matrix temp;
      temp.Product(tmp,R);
      for(int i=0;i<n;i++)
      {
	 for(int j=0;j<n;j++)
	    tmp[i][j] = temp[i][j] * sqrt(weight[j]);
      }
      final.Product(tmp,tmp2);
      svd.Decompose(final);
      for(int i=0;i<n;i++)
	 lambda[i] = fabs(svd.w[i]);
   }
}

void Meta::BurdenAssoc(String method,GroupFromAnnotation & group,Vector *& maf,Vector *& stats,Vector *& cond_stats,Matrix *& cov,Matrix *& cond_cov,Vector *& singleEffSize,Vector *& singlePvalue)
{
   printf("Performing %s tests ...\n",method.c_str());
   //calculate final results here

   Vector pvalue_burden,pvalue_burden_cond;

   IFILE output;
   String filename;
   if(prefix =="")
      filename = "meta."+method + ".results";
   else if(prefix.Last()=='.' || prefix.Last()=='/')
      filename = prefix +  "meta."+method +".results";
   else
      filename = prefix + ".meta."+method +".results";
   output=ifopen(filename,"w",InputFile::UNCOMPRESSED);

   String method_out = method;
   method_out +="_";
   method_out += MAF_cutoff;

   IFILE reportOutput;
   String reportFile;
   if(report)
   {
      if(prefix =="")
	 reportFile = "meta.tophits."+method+".tbl";
      else if(prefix.Last()=='.' || prefix.Last()=='/')
	 reportFile = prefix +  "meta.tophits."+method+".tbl";
      else
	 reportFile = prefix + ".meta.tophits."+method+".tbl";
      reportOutput=ifopen(reportFile,"w",InputFile::UNCOMPRESSED);
      if(cond!="")
	 ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\tCOND_EFFSIZES\tCOND_PVALUES\n");
      else
	 ifprintf(reportOutput,"GENE\tMETHOD\tGENE_PVALUE\tMAF_CUTOFF\tACTUAL_CUTOFF\tVARS\tMAFS\tEFFSIZES\tPVALUES\n");
   }
   ifprintf(output,"##Method=Burden\n");
   ifprintf(output,"##STUDY_NUM=%d\n",studies.Length());
   ifprintf(output,"##TotalSampleSize=%d\n",total_N);
   if(cond!="")
      if(fullResult)
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
      else
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\tCOND_EFFSIZE\tCOND_PVALUE\n");
   else
      if(fullResult)
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
      else
	 ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tPVALUE\n");
   double numerator,denominator,chisq,pvalue,cond_num,cond_denom=_NAN_;
   StringArray chr_plot,geneLabels;
   Vector pos_plot;
   for(int g=0;g<group.annoGroups.Length();g++)
   {
      if(maf[g].Length()==0)
      {
	 if(cond!="")
	    if(fullResult)
	       ifprintf(output,"%s\t0\t-\t-\t-\t-\t0\t0\t0\t-\t-\t-\t-\n",group.annoGroups[g].c_str());
	    else
	       ifprintf(output,"%s\t0\t-\t0\t0\t0\t-\t-\t-\t-\n",group.annoGroups[g].c_str());
	 else
	    if(fullResult)
	       ifprintf(output,"%s\t0\t-\t-\t-\t-\t0\t0\t0\t-\t-\n",group.annoGroups[g].c_str());
	    else
	       ifprintf(output,"%s\t0\t-\t0\t0\t0\t-\t-\n",group.annoGroups[g].c_str());
	 continue;
      }
      double average_af = maf[g].Average();
      double min_af = maf[g].Min();
      double max_af = maf[g].Max();

      String var;
      for(int i=0;i<maf[g].Length()-1;i++)
      {
	 var += group.SNPlist[g][i] + ";";
      }
      var += group.SNPlist[g][maf[g].Length()-1];

      Vector weight;
      weight.Dimension(maf[g].Length());

      if(method=="burden")
	 for(int w=0;w<weight.Length();w++)
	    weight[w] = 1.0;
      else if(method=="MB")
	 for(int w=0;w<weight.Length();w++)
	    weight[w] = 1.0/sqrt(maf[g][w]*(1.0-maf[g][w]));

      numerator  = weight.InnerProduct(stats[g]);
      Vector tmp;
      tmp.Dimension(group.SNPlist[g].Length());

      for(int i=0;i<tmp.Length();i++)
      {
	 tmp[i] = weight.InnerProduct(cov[g][i]);
      }
      denominator = tmp.InnerProduct(weight);

      //printf("num %g, denom %g\n",numerator,denominator);

      if(cond!="")
      {
	 cond_num = weight.InnerProduct(cond_stats[g]);
	 for(int i=0;i<tmp.Length();i++)
	 {
	    tmp[i] = weight.InnerProduct(cond_cov[g][i]);
	 }
	 cond_denom = tmp.InnerProduct(weight);
      }
      /*
	 if(group.annoGroups[g]=="KIAA1024")
	 {
	 printf("cond_num %g, cond_denom %g\n",cond_num,cond_denom);
	 printf("cond_stats are:\n");
	 for(int i=0;i<cond_stats[g].Length();i++)
	 printf("%g\t",cond_stats[g][i]);
	 printf("\n");
	 printf("stats are:\n");
	 for(int i=0;i<stats[g].Length();i++)
	 printf("%g\t",stats[g][i]);
	 printf("\n");
	 printf("covs are:\n");
	 for(int i=0;i<cov[g].rows;i++)
	 {
	 for(int j=0;j<cov[g].cols;j++)
	 printf("%g\t",cov[g][i][j]);
	 printf("\n");
	 }
	 printf("cond covs are:\n");
	 for(int i=0;i<cond_cov[g].rows;i++)
	 {
	 for(int j=0;j<cond_cov[g].cols;j++)
	 printf("%g\t",cond_cov[g][i][j]);
	 printf("\n");
	 }
	 }
       */
      if(denominator==0.0)
      {
	 if(fullResult)
	 {
	    ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());

	    for(int i=0;i<maf[g].Length()-1;i++)
	       ifprintf(output,"%g,",maf[g][i]);
	    ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);

	    for(int i=0;i<singleEffSize[g].Length()-1;i++)
	       ifprintf(output,"%g,",singleEffSize[g][i]);
	    ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);

	    for(int i=0;i<singlePvalue[g].Length()-1;i++)
	       ifprintf(output,"%g,",singlePvalue[g][i]);
	    ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);
	    if(cond!="")
	       ifprintf(output,"%g\t%g\t%g\t-\t-\t-\t-\n",average_af,min_af,max_af);
	    else
	       ifprintf(output,"%g\t%g\t%g\t-\t-\n",average_af,min_af,max_af);
	 }
	 else
	    if(cond!="")
	       ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t-\t-\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af);
	    else
	       ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t-\t-\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af);
	 continue;
      }

      chisq = numerator*numerator/denominator;
      pvalue = pchisq(chisq,1,0,0);
      double effSize = numerator/denominator;
      double cond_chisq,cond_effSize,cond_pvalue;
      /*
	 if(group.annoGroups[g]=="KIAA1024")
	 printf("gene:%s, numerator: %g,denominator %g,chisq %g\n",group.annoGroups[g].c_str(),numerator,denominator,chisq);
	 printf("numerators:\n");
	 for(int i=0;i<stats[g].dim;i++)
	 printf("stats are: %g ",stats[g][i]);
	 printf("numerator is:%g\n",numerator);

	 printf("denominators:\n");
	 for(int i=0;i<cov[g].rows;i++)
	 {
	 for(int j=0;j<cov[g].cols;j++)
	 printf("covs are: %g ",cov[g][i][j]);
	 printf("\n");
	 }

	 printf("denominators are:%g,pvalue is: %g\n",denominator,pvalue);
       */
      bool disect=false;
      while(pvalue==0.0)
      {
	 disect=true;
	 chisq *= 0.999;
	 pvalue = pchisq(chisq,1,0,0);
      }
      bool cond_disect =false;
      if(cond!="")
      {
	 if(cond_denom==0)
	 {
	    cond_effSize =0.0;
	    cond_pvalue =1.0;;
	 }
	 else
	 {
	    cond_chisq = cond_num*cond_num/cond_denom;
	    cond_effSize = cond_num/cond_denom;
	    cond_pvalue = pchisq(cond_chisq,1,0,0);
	    while(cond_pvalue==0.0)
	    {
	       cond_disect=true;
	       cond_chisq *= 0.999;
	       cond_pvalue = pchisq(cond_chisq,1,0,0);
	    }
	 }
	 pvalue_burden_cond.Push(cond_pvalue);
      }

      if(fullResult)
      {
	 ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str());

	 for(int i=0;i<maf[g].Length()-1;i++)
	    ifprintf(output,"%g,",maf[g][i]);
	 ifprintf(output,"%g\t",maf[g][maf[g].Length()-1]);

	 for(int i=0;i<singleEffSize[g].Length()-1;i++)
	    ifprintf(output,"%g,",singleEffSize[g][i]);
	 ifprintf(output,"%g\t",singleEffSize[g][singleEffSize[g].Length()-1]);

	 for(int i=0;i<singlePvalue[g].Length()-1;i++)
	    ifprintf(output,"%g,",singlePvalue[g][i]);
	 ifprintf(output,"%g\t",singlePvalue[g][singlePvalue[g].Length()-1]);

	 if(cond!="")
	    ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
	 else
	    ifprintf(output,"%g\t%g\t%g\t%g\t%s%g\n",average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
      }
      else
      {
	 if(cond!="")
	    ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue,cond_effSize,cond_disect?"<":"",cond_pvalue);
	 else
	    ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%s%g\n",group.annoGroups[g].c_str(),group.SNPlist[g].Length(),var.c_str(),average_af,min_af,max_af,effSize,disect?"<":"",pvalue);
      }
      if(pvalue <report_pvalue_cutoff && report)
      {
	 StringArray variants;
	 variants.AddTokens(var,";");
	 if(cond!="")
	    for(int v=0;v<maf[g].Length();v++)
	       ifprintf(reportOutput,"%s\t%s\t%s%g\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,cond_disect?"<":"",cond_pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
	 else
	    for(int v=0;v<maf[g].Length();v++)
	       ifprintf(reportOutput,"%s\t%s\t%s%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),disect?"<":"",pvalue,MAF_cutoff,MAF_cutoff,variants[v].c_str(),maf[g][v],singleEffSize[g][v],singlePvalue[g][v]);
      }
      pvalue_burden.Push(pvalue);
      geneLabels.Push(group.annoGroups[g]);
      StringArray tmp_SNPname;
      tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
      chr_plot.Push(tmp_SNPname[0]);
      pos_plot.Push(tmp_SNPname[1].AsDouble());
   }
   String name = method;
   name += " (maf<";
   name += MAF_cutoff;
   name +=  ")";
   String extraname = "";
   String demo;
   if(pvalue_burden.Length()>0)
   {
      writepdf.Draw(pdf,geneLabels,pvalue_burden,chr_plot,pos_plot,name,extraname,demo,true);
      if(cond!="")
      {
	 name += "conditional analysis";
	 writepdf.Draw(pdf,geneLabels,pvalue_burden_cond,chr_plot,pos_plot,name,extraname,demo,true);
      }
   }
   ifclose(output);
   if(report)
      ifclose(reportOutput);
   printf("  done.\n\n");
}

double Meta::GetBetaDensity(double a, double b, double x)
{
   double density;
   density = exp(gammln(a+b)-gammln(a)-gammln(b)+(a-1.0)*log(x)+(b-1.0)*log(1.0-x));
   return density;
}

double Meta::CalculateCorrCoef(Vector & a,Vector & b)
{
   double coef;
   int n = a.Length();
   double sum_a = a.Sum();
   double sum_b = b.Sum();
   coef =  n*a.InnerProduct(b) - sum_a*sum_b;
   coef /= sqrt((n*a.SumSquares()-sum_a*sum_a)*(n*b.SumSquares()-sum_b*sum_b));  
   return coef;
}

//calculate the genotype covariance between G and X
void Meta::CalculateGenotypeCov(SummaryFileReader & covReader,String chr, int pos,int study,Vector & result)
{
   result.Dimension(commonVar_study[study].Length(),0.0);
   for(int cv=0;cv<commonVar_study[study].Length();cv++)
   {
      if (chr != common_chr[commonVar_study[study][cv]])
	 continue;
      if (abs(common_pos[commonVar_study[study][cv]] - pos) > 1000000)
	 continue;
      if(common_pos[commonVar_study[study][cv]] > pos)
      {
	 if(!covReader.ReadRecord(chr,pos))
	    continue;
	 StringArray tmp;
	 tmp.AddTokens(covReader.marker_nearby,",");
	 String pos_str;
	 pos_str = common_pos[commonVar_study[study][cv]];
	 int loc = tmp.Find(pos_str);
	 if(loc==-1)
	    continue;
	 tmp.Clear();
	 tmp.AddTokens(covReader.marker_cov,",");
	 result[cv] =  tmp[loc].AsDouble();
      }
      else
      {
	 int loc = commonVar_markers_in_window[study][cv].Find(pos);
	 if(loc==-1)
	    continue;
	 result[cv] = commonVar_marker_cov[study][cv][loc];
	 /*
	    if(pos==79760530)
	    printf("chr %s pos %d GXcov %g\n",chr.c_str(),pos,commonVar_marker_cov[study][cv][loc]);
	  */    }
   }
}

double Meta::GrabGenotypeCov(SummaryFileReader & covReader,int study,String chr1,String pos1,String chr2,String pos2,String & SNP1, String & SNP2)
{
   double result = 0.0;
   if(chr1!=chr2)
      return result;

   int skip = SNPexclude.Integer(study+":"+chr1+":"+pos1);
   if(skip!=-1)
      return result;
   skip = SNPexclude.Integer(study+":"+chr2+":"+pos2);
   if(skip!=-1)
      return result;

   StringArray flip,tmp;
   int marker_idx;

   if(pos1.AsInteger()<pos2.AsInteger())
   {
      //Check in current hash table first
      marker_idx = covReader.markerPosHash.Integer(chr1 + ":"+pos1);
      //if this record has been hashed
      if(marker_idx!=-1)
      {
	 tmp.AddTokens(covReader.markerNearby[marker_idx],",");
	 int loc = tmp.Find(pos2);
	 if(loc==-1)
	    return result;
	 tmp.Clear();
	 tmp.AddTokens(covReader.markerNearbyCov[marker_idx],",");
	 result = tmp[loc].AsDouble();
	 return result;
      }
      if(!covReader.ReadRecord(chr1,pos1.AsInteger()))
	 return result;
      tmp.Clear();
      tmp.AddTokens(covReader.marker_nearby,",");
      int loc = tmp.Find(pos2);
      if(loc==-1)
	 return result;
      tmp.Clear();
      tmp.AddTokens(covReader.marker_cov,",");
      result = tmp[loc].AsDouble();
   }
   else
   {
      //Check in current hash table first
      marker_idx = covReader.markerPosHash.Integer(chr2 + ":"+pos2);
      //if this record has been hashed
      if(marker_idx!=-1)
      {
	 StringArray tmp;
	 tmp.AddTokens(covReader.markerNearby[marker_idx],",");
	 int loc = tmp.Find(pos1);
	 if(loc==-1)
	    return result;
	 tmp.Clear();
	 tmp.AddTokens(covReader.markerNearbyCov[marker_idx],",");
	 result = tmp[loc].AsDouble();
	 return result;
      }
      if(!covReader.ReadRecord(chr2,pos2.AsInteger()))
	 return result;
      tmp.Clear();
      tmp.AddTokens(covReader.marker_nearby,",");
      int loc = tmp.Find(pos1);
      if(loc==-1)
	 return result;
      tmp.Clear();
      tmp.AddTokens(covReader.marker_cov,",");
      result = tmp[loc].AsDouble();
   }

   //if this marker is flipped then markers from the entire row
   //should have cov multiply by -1.0.
   double factor1=1.0,factor2=1.0;
   if(flipSNP.Integer(study+":"+chr1+":"+pos1)!=-1)
      factor1=-1.0;
   if(flipSNP.Integer(study+":"+chr2+":"+pos2)!=-1)
      factor2=-1.0;
   result *= (factor1*factor2);

   return result;
}

void Meta::CalculateXXCov(int study,Matrix & result)
{
   for(int i=0;i<result.rows;i++)
   {
      result[i][i] = commonVar_marker_cov[study][i][0];
      for(int j=i+1;j<result.cols;j++)
      {
	 if(common_chr[commonVar_study[study][j]] != common_chr[commonVar_study[study][i]])
	    continue;
	 if(common_pos[commonVar_study[study][i]]<common_pos[commonVar_study[study][j]])
	 {
	    int loc = commonVar_markers_in_window[study][i].Find(common_pos[commonVar_study[study][j]]);
	    if(loc==-1)
	       continue;
	    result[i][j] = result[j][i] = commonVar_marker_cov[study][i][loc];
	 }
	 else
	 {
	    int loc = commonVar_markers_in_window[study][j].Find(common_pos[commonVar_study[study][i]]);
	    if(loc==-1)
	       continue;
	    result[i][j] = result[j][i] = commonVar_marker_cov[study][j][loc];
	 }
      }
   }
   SVD svd;
   svd.InvertInPlace(result);
}

void Meta::RevertAllele(String SNP, String & newSNP)
{
   StringArray tmp;
   tmp.AddTokens(SNP,":");
   newSNP = tmp[0]+":"+tmp[1]+":"+tmp[3]+":"+tmp[2];
}

