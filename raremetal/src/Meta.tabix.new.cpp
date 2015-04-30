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

String Meta::studyName = "";
String Meta::prefix = ""; 
bool Meta::Burden = false;
bool Meta::MB = false;
bool Meta::SKAT = false;
bool Meta::VTa = false;
bool Meta::VTp = false;
bool Meta::outvcf =false; 
bool Meta::fullResult =false; 
bool Meta::report = false; 
double  Meta::report_pvalue_cutoff = 1e-06; 
bool Meta::founderAF = false; 
double Meta::HWE = 0.0;
double Meta::CALLRATE = 0.0;
double Meta::MAF_cutoff = 0.05;
int Meta::marker_col = 2;
int Meta::cov_col = 3;
String Meta::cond = "";

Meta::Meta(){}
Meta::~Meta()
{}

//This function read all study names from a file 
//and save the names in the StringArray files
void Meta::Prepare()
{
   String filename;
   if(prefix != "")
   {
      if(prefix.FastFindLastChar('.')!=-1 || prefix.FastFindLastChar('/')!=-1)
	 filename = prefix + ".meta.plots.pdf";
      else
	 filename = prefix + ".meta.plots.pdf";
   }
   else
      filename = "meta.plots.pdf";
   pdf.OpenFile(filename);

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
	 if(buffer.FindChar('#')!=-1)
	    continue;
	 commonVar.AddTokens(buffer, " \t\n");
	 for(int i=0;i<commonVar.Length();i++)
	 {
	    StringArray tmp;
	    tmp.AddTokens(commonVar[i],":");
	    common_chr.Push(tmp[0]);
	    common_pos.Push(tmp[1]);
	    common_ref.Push(tmp[2]);
	    common_alt.Push(tmp[3]);
	 }
      }
      ifclose(condFile);

      int numStudy = studies.Length();
      commonVar_study = new IntArray [numStudy];
      commonVar_effSize = new Vector [numStudy];
      commonVar_U = new Vector [numStudy];
      commonVar_V = new Vector [numStudy];
      commonVar_markers_in_window = new IntArray * [numStudy];
      commonVar_marker_cov = new Vector * [numStudy];
      XX_inv = new Matrix [numStudy];
      cond_status = new bool [numStudy];

      for(int s=0;s<studies.Length();s++)
      {
	 SummaryFileReader statReader,covReader;

	 String filename = studies[s] + ".singlevar.score.txt.gz";
	 String covFilename = studies[s] + ".singlevar.cov.txt.gz";

	 if(!statReader.ReadTabix(filename))
	 {
	    filename = studies[s] + ".MetaScore.assoc.gz";
	    if(!statReader.ReadTabix(filename))
	    {
	       error("Can not file file %s. Please check file name for study %s\n",filename.c_str(),studies[s].c_str());
	    }
	 }
	 if(!covReader.ReadTabix(covFilename))
	 {
	    covFilename = studies[s] + ".MetaCov.assoc.gz";
	    if(!covReader.ReadTabix(covFilename))
	    {
	       error("Can not open file %s. Please check file name for study %s\n",covFilename.c_str(),studies[s].c_str());
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
	       commonVar_study[s].Push(i);
	       commonVar_effSize[s].Push(record[15].AsDouble());
	       commonVar_U[s].Push(record[13].AsDouble());
	       double v = record[14].AsDouble();
	       commonVar_V[s].Push(v*v);
	       commonVar_betaHat = new Vector [numStudy];
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
	       error("Covariance file is not consistent with summary statistics file.\n");
	    StringArray tmp,tmp_markers;
	    tmp.AddTokens(covReader.buffer,"\t");
	    tmp_markers.AddTokens(tmp[2],",");
	    for(int j=0;j<tmp_markers.Length();j++)
	       commonVar_markers_in_window[s][i].Push(tmp_markers[j].AsInteger());
	    tmp_markers.Clear();
	    tmp_markers.AddTokens(tmp[3],",");
	    for(int j=0;j<tmp_markers.Length();j++)
	       commonVar_marker_cov[s][i].Push(tmp_markers[j].AsDouble());
	 }
	 CalculateXXCov(s,XX_inv[s]);
	 //calculate beta hat for common variants
	 commonVar_betaHat[s].Dimension(dim);
	 for(int i=0;i<commonVar_betaHat[s].dim;i++)
	    commonVar_betaHat[s][i] = XX_inv[s][i].InnerProduct(commonVar_U[s]);
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
void Meta::PoolSummaryStat()
{
   //usefulSize and usefulAC have the pooled N and AC information 
   StringIntHash usefulSize;
   StringDoubleHash usefulAC; 
   //refalt use study:chr:pos as key and ref:alt as value
   StringHash refalt;
   StringIntHash SNP_DirectionByStudy,SNP_DirectionByStudyNoAllele;
   StringArray directions;

   total_N=0;
   int flipCount=0;

   for(int study=0;study<studies.Length();study++)
   {
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
	 covReader.counter=0;
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
	 error("Cannot open file: %s!\n",filename.c_str());

      String buffer;
      StringArray tokens;
      int current_N =0;
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
		  for(int s=0;s<study+1-len;s++)
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
		  for(int s=0;s<study+1-len;s++)
		     directions[directions_idx] += "?";
	       }
	       directions[directions_idx] += "!";
	    }
	    else
	    {
	       String * refalt_saved = (String *) refalt.Object(tokens[0]+":"+tokens[1]);
	       if(refalt_saved != NULL)
	       {
		  //if this variant at this chr:pos has been hashed but with different ref/alt alleles
		  //then skip this variant
		  SNPexclude.SetInteger(study+":"+tokens[0]+":"+tokens[1],0);
		  directions_idx = SNP_DirectionByStudyNoAllele.Integer(tokens[0]+":"+tokens[1]);
		  if(directions[directions_idx].Length()<study)
		  {
		     int len = directions[directions_idx].Length();
		     for(int s=0;s<study+1-len;s++)
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

	    //update SNP_DirectionByStudy
	    int directions_idx;
	    String direction = "+";
	    if(tokens[15-adjust]=="NA" || tokens[15-adjust]=="na" || tokens[15-adjust]=="nan" || tokens[15-adjust]=="-")
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
	       for(int s=0;s<study+1-len;s++)
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
	       SNPexclude.SetInteger(study+":"+tokens[0]+":"+tokens[1],0);
	       int directions_idx = SNP_DirectionByStudyNoAllele.Integer(tokens[0]+":"+tokens[1]);
	       if(directions[directions_idx].Length()<study)
	       {
		  int len = directions[directions_idx].Length();
		  for(int s=0;s<study+1-len;s++)
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

	    if(tokens[15-adjust]=="NA" || tokens[15-adjust]=="na" || tokens[15-adjust]=="nan" || tokens[15-adjust]=="-")
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
   }

   //calculate pooled allele frequencies
   for(int i=0;i<usefulAC.Capacity();i++)
   {
      if(!usefulAC.SlotInUse(i)) 
      {
	 continue;
      }
      String SNPname = usefulAC[i];
      double AC = usefulAC.Double(i);
      int N = usefulSize.Integer(SNPname);
      double maf;
      if(founderAF)
      {
	 maf = AC/(2.0*total_N);
      }
      else
      {
	 maf = AC/(2.0*N);
      }
      if(maf>0.5) 
	 maf = 1.0-maf;

      SNPmaf.SetDouble(SNPname,maf);
   } //now all SNP maf are hashed in SNPmaf.

   //calculate single variant meta results
   printf("Performing Single variant meta analysis ...\n");
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

   IFILE vcfout;
   if(outvcf)
   {
      String filename;
      if(prefix=="")
	 filename = "pooled.variants.vcf";
      else if(prefix.Last()=='.' || prefix.Last()=='/')
	 filename = prefix + "pooled.variants.vcf";
      else 
	 filename = prefix + ".pooled.variants.vcf";

      vcfout = ifopen(filename,"w",InputFile::UNCOMPRESSED);
      ifprintf(vcfout,"#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
   }

   Vector pvalueAll,pvalue1,pvalue5,pvalueAll_cond,pvalue1_cond,pvalue5_cond;
   StringArray chr_plot;
   Vector pos_plot;

   //TODO: need to output variants sorted by chr and pos
   for(int i=0;i<SNPmaf.Capacity();i++)
   {
      if(!SNPmaf.SlotInUse(i))
      {
	 continue;
      }
      String SNPname = SNPmaf[i];
      double U = SNPstat.Double(SNPname);
      double V = SNP_Vstat.Double(SNPname);
      double maf = SNPmaf.Double(SNPname);
      int direction_idx = SNP_DirectionByStudy.Integer(SNPname);
      if(directions[direction_idx].Length()<studies.Length())
      {
	 int len = directions[direction_idx].Length();
	 for(int s=0;s<studies.Length()-len;s++)
	    directions[direction_idx] += "?";
      }

      String direction = directions[direction_idx];
      //printf("direction is %s,SNPname is %s\n",directions[direction_idx].c_str(), SNPname.c_str());

      if(maf>0 && maf<1.0 && V!=0.0)
      {
	 double chisq = U*U/V;
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
	    pvalue1.Push(pvalue);
	 if(maf<0.05)
	    pvalue5.Push(pvalue);
	 chr_plot.Push(tokens[0]);
	 pos_plot.Push(tokens[1].AsInteger());

	 if(cond!="")
	 {
	    double cond_U = SNPstat_cond.Double(SNPname);
	    double cond_V = SNP_Vstat_cond.Double(SNPname);
	    double chisq = cond_U*cond_U/cond_V;;
	    cond_pvalue = pchisq(chisq,1,0,0);
	    cond_effSize = cond_U/cond_V;
	    bool cond_disect=false;
	    while(cond_pvalue==0.0)
	    {
	       disect=true;
	       chisq *= 0.999;
	       cond_pvalue = pchisq(chisq,1,0,0);
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
      }
   }

   writepdf.Draw(pdf,pvalueAll,pvalue1,pvalue5,chr_plot,pos_plot,"single variant results");
   if(cond!="")
      writepdf.Draw(pdf,pvalueAll_cond,pvalue1_cond,pvalue5_cond,chr_plot,pos_plot,"conditonal analysis");

   printf("  done.\n\n");

   ifclose(output);
   if(outvcf) ifclose(vcfout);
}

void Meta::Run(GroupFromAnnotation & group)
{
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
	 double af = SNPmaf.Double(group.SNPlist[g][m]);
	 double singleP = singleVarPvalue.Double(group.SNPlist[g][m]);
	 double singleEff = singleVarEff.Double(group.SNPlist[g][m]);
	 if(af==_NAN_)
	 {
	    flipStatus=true;
	    StringArray SNPflip;
	    SNPflip.AddTokens(group.SNPlist[g][m],":");
	    newSNP = SNPflip[0]+":"+SNPflip[1]+":"+SNPflip[3]+":"+SNPflip[2];
	    af = SNPmaf.Double(newSNP);
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
      //printf("there are %d SNPs in this group\n",count);
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
      String covFilename = studies[study] + ".singlevar.cov.txt.gz";
      marker_col = 2;
      cov_col = 3;

      if(!covReader.ReadTabix(covFilename))
      {
	 marker_col = 4;
	 cov_col = 5;
	 covFilename = studies[study] + ".MetaCov.assoc.gz";
	 if(!covReader.ReadTabix(covFilename))
	 {
	    error("Can not open file %s. Please check file name for study %s\n",covFilename.c_str(),studies[study].c_str());
	 }
      }

      //   printf("Updating group stats ...\n");
      //update group statistics
      for(int g=0;g<group.annoGroups.Length();g++)
      {
	 printf("doing group %d\n",g);
	 int count = group.SNPlist[g].Length();
	 StringArray chr,pos;
	 for(int i=0;i<count;i++)
	 {
	    StringArray tmp;
	    tmp.AddTokens(group.SNPNoAllele[g][i],":");
	    chr.Push(tmp[0]);
	    pos.Push(tmp[1]);
	 } //now pos has all the positions of markers in group g.

	 Matrix cov_i,GX,extra_cov_i;
	 cov_i.Dimension(count,count,0.0);
	 if(cond!="")
	 {
	    extra_cov_i.Dimension(count,count,0.0);
	    GX.Dimension(count,XX_inv[study].cols);
	    //GX_beta.Dimension(count);
	 }
	 for(int m=0;m<count;m++)
	 {
	    if(cond!="")
	    {
	       CalculateGenotypeCov(covReader,chr[m],pos[m].AsInteger(),study,GX[m]);
	       //GX_beta[m] = GX[m].InnerProduct(commonVar_betaHat[study]);
	       //cond_stats[g][m] = stats[g][m]-GX_beta[m];
	       for(int c=0;c<count;c++)
		  cond_stats[g][m][c] = SNPstat_cond.Double(group.SNPlist[g][c]);
	    }
	    for(int s=m;s<count;s++)
	    {
	       cov_i[m][s] = GrabGenotypeCov(covReader,study,chr[m],pos[m],chr[s],pos[s],group.SNPlist[g][m],group.SNPlist[g][s]);
	    }
	 }
	 cov[g].Add(cov_i);
	 if(cond!="")
	 {
	    Matrix GX_trans,tmp;
	    GX_trans.Transpose(GX);
	    tmp.Product(GX,XX_inv[study]);
	    extra_cov_i.Product(tmp,GX_trans);
	    extra_cov_i.Multiply(-1.0);
	    cond_cov[g].Add(extra_cov_i);
	 }
      }
   }

   for(int g=0;g<group.annoGroups.Length();g++)
   {
      for(int r=0;r<cov[g].rows;r++)
	 for(int c=r+1;c<cov[g].cols;c++)
	 {
	    cov[g][c][r] = cov[g][r][c];
	    if(cond!="")
	       cond_cov[g][c][r] = cond_cov[g][r][c];
	 }
      if(cond!="")
	 cond_cov[g].Add(cov[g]);
   }

   String method = "";
   if(Burden)
   {
      method = "Burden";
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
      StringArray chr_plot;

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
	 }

	 pvalue_VT.Push(pvalue);

	 StringArray tmp_SNPname;
	 tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
	 chr_plot.Push(tmp_SNPname[0]);
	 pos_plot.Push(tmp_SNPname[1].AsDouble());

	 if(lower) delete [] lower;
	 if(upper) delete [] upper;
	 if(mean) delete [] mean;
      }

      String name = "VT";
      String extraname = "";
      extraname += MAF_cutoff;
      writepdf.Draw(pdf,pvalue_VT,chr_plot,pos_plot,name,extraname);


      ifclose(output);
      if(report)
	 ifclose(reportOutput);
      printf("Finished Variable Threshold tests.\n\n");
   }

   if(SKAT)
   {
      printf("Performing SKAT ...\n");
      //calculate Q statistics here
      Vector pvalue_SKAT,pos_plot,cond_pvalue_SKAT;
      StringArray chr_plot;
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
      SVD svd;
      Matrix Pmix;
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
	 double cond_Qstat;
	 if(cond!="")
	    cond_Qstat = cond_tmp.InnerProduct(cond_stats[g]);
	 //printf("Qstat is: %g\n",Qstat);
	 svd.Decompose(cov[g]);
	 Matrix L;
	 L.Dimension(n,n);
	 Matrix Pmix_tmp,Pmix;
	 Pmix.Dimension(n,n);
	 Pmix_tmp.Dimension(n,n);
	 Matrix L_trans;

	 for(int j=0;j<n;j++)
	 {
	    double sqrt_w = sqrt(fabs(svd.w[j]));
	    for(int i=0;i<n;i++)
	    {
	       L[i][j] = svd.u[i][j] * sqrt_w;
	    }
	 }
	 for(int j=0;j<n;j++)
	 {
	    for(int i=0;i<n;i++)
	    {
	       Pmix_tmp[i][j] = L[i][j] * weight[j];
	    }
	 }
	 L_trans.Transpose(L);
	 Pmix.Product(Pmix_tmp,L_trans);
	 svd.Decompose(Pmix);
	 double * lambda = new double [n];
	 for(int i=0;i<n;i++)
	 {
	    lambda[i] = fabs(svd.w[i]);
	 }

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
	 disect_itr=0;
	 while((pvalue<=0.0 ||pvalue==2.0 || isnan(pvalue))  && disect_itr<10000)
	 {
	    Qstat_dav*=0.99;
	    pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
	    disect_itr++;
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
	 disect_itr=0;
	 while((pvalue_liu<=0.0 ||pvalue_liu==2.0 || isnan(pvalue_liu))  && disect_itr<10000)
	 {
	    Qstat_liu*=0.99;
	    pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	    disect_itr++;
	 }

	 if(cond!="")
	 {
	    svd.Decompose(cond_cov[g]);
	    for(int j=0;j<n;j++)
	    {
	       double sqrt_w = sqrt(fabs(svd.w[j]));
	       for(int i=0;i<n;i++)
	       {
		  L[i][j] = svd.u[i][j] * sqrt_w;
	       }
	    }
	    for(int j=0;j<n;j++)
	    {
	       for(int i=0;i<n;i++)
	       {
		  Pmix_tmp[i][j] = L[i][j] * weight[j];
	       }
	    }
	    L_trans.Transpose(L);
	    Pmix.Product(Pmix_tmp,L_trans);
	    svd.Decompose(Pmix);
	    for(int i=0;i<n;i++)
	    {
	       lambda[i] = fabs(svd.w[i]);
	    }

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
	    disect_itr=0;
	    while((cond_pvalue<=0.0 ||cond_pvalue==2.0 || isnan(cond_pvalue))  && disect_itr<10000)
	    {
	       Qstat_dav*=0.99;
	       cond_pvalue = MixChidist(lambda, n, Qstat_dav,"Davies");
	       disect_itr++;
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
	    disect_itr=0;
	    while((cond_pvalue_liu<=0.0 ||cond_pvalue_liu==2.0 || isnan(cond_pvalue_liu))  && disect_itr<10000)
	    {
	       Qstat_liu*=0.99;
	       cond_pvalue_liu = MixChidist(lambda, n, Qstat_liu,"Liu");
	       disect_itr++;
	    }
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
	 if((pvalue <report_pvalue_cutoff || pvalue_liu<report_pvalue_cutoff)  && report)
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
      }

      String name = "SKAT";
      String extraname = "";
      extraname += MAF_cutoff;
      writepdf.Draw(pdf,pvalue_SKAT,chr_plot,pos_plot,name,extraname);
      if(cond!="")
	 writepdf.Draw(pdf,cond_pvalue_SKAT,chr_plot,pos_plot,name,extraname);

      ifclose(output);
      if(report)
	 ifclose(reportOutput);
      printf("Finished working on SKAT.\n\n");
   }

   pdf.CloseFile();

   //housekeeping
   if(maf) delete [] maf;
   if(stats) delete [] stats;
   if(cov) delete [] cov;
   if(singlePvalue) delete [] singlePvalue;
   if(singleEffSize) delete [] singleEffSize;
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
   StringArray chr_plot;
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

      if(method=="Burden")
	 for(int w=0;w<weight.Length();w++)
	    weight[w] = 1.0;
      else if(method=="MB")
	 for(int w=0;w<weight.Length();w++)
	    weight[w] = 1.0;

      numerator  = weight.InnerProduct(stats[g]);
      Vector tmp;
      tmp.Dimension(group.SNPlist[g].Length());

      for(int i=0;i<tmp.Length();i++)
      {
	 tmp[i] = weight.InnerProduct(cov[g][i]);
      }
      denominator = tmp.InnerProduct(weight);
      printf("num %g, denom %g\n",numerator,denominator);

      if(cond!="")
      {
	 cond_num = weight.InnerProduct(cond_stats[g]);
	 tmp.Zero();
	 for(int i=0;i<tmp.Length();i++)
	 {
	    tmp[i] = weight.InnerProduct(cond_cov[g][i]);
	 }
	 cond_denom = tmp.InnerProduct(weight);
      }

      if(denominator==0.0 || cond_denom==0.0)
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
      //printf("gene:%s, numerator: %g,denominator %g,chisq %g\n",group.annoGroups[g].c_str(),numerator,denominator,chisq);

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
	 cond_chisq = cond_num*cond_num/cond_denom;
	 cond_effSize = cond_num/cond_denom;
	 cond_pvalue = pchisq(cond_chisq,1,0,0);
	 while(cond_pvalue==0.0)
	 {
	    cond_disect=true;
	    cond_chisq *= 0.999;
	    cond_pvalue = pchisq(cond_chisq,1,0,0);
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
      StringArray tmp_SNPname;
      tmp_SNPname.AddTokens(group.SNPlist[g][0],":");
      chr_plot.Push(tmp_SNPname[0]);
      pos_plot.Push(tmp_SNPname[1].AsDouble());
   }
   String name = method;
   String extraname = "";
   extraname += MAF_cutoff;
   if(pvalue_burden.Length()>0)
   {
      writepdf.Draw(pdf,pvalue_burden,chr_plot,pos_plot,name,extraname);
      if(cond!="")
	 writepdf.Draw(pdf,pvalue_burden_cond,chr_plot,pos_plot,name,extraname);
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
   result.Dimension(commonVar_study[study].Length());
   result.Zero();

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
	 StringArray currentRecord,tmp;
	 currentRecord.AddTokens(covReader.buffer," \t");
	 tmp.AddTokens(currentRecord[2],",");
	 String pos_str = common_pos[commonVar_study[study][cv]];
	 int loc = tmp.Find(pos_str);
	 if(loc==-1)
	    continue;
	 tmp.Clear();
	 tmp.AddTokens(currentRecord[3],",");
	 result[cv] =  tmp[loc].AsDouble()*SampleSize[study];
      }
      else
      {
	 int loc = commonVar_markers_in_window[study][cv].Find(pos);
	 if(loc==-1)
	    continue;
	 result[cv] = commonVar_marker_cov[study][cv][loc]*SampleSize[study];
      }
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

   StringArray flip;
   int marker_idx;

   if(pos1.AsInteger()<pos2.AsInteger())
   {
      //Check in current hash table first
      marker_idx = covReader.markerPosHash.Integer(chr1 + ":"+pos1);
      //if this record has been hashed
      if(marker_idx!=-1)
      {
	 StringArray tmp;
	 tmp.AddTokens(covReader.markerNearby[marker_idx],",");
	 int loc = tmp.Find(pos);
	 if(loc==-1)
	    return result;
	 tmp.Clear();
	 tmp.AddTokens(covReader.markerNearbyCov[marker_idx],",");
	 result = tmp[loc].AsDouble()*SampleSize(study);
	 return result;
      }
      if(!covReader.ReadRecord(chr1,pos1.AsInteger()))
	 return result;
      StringArray tmp;
      tmp.AddTokens(covReader.marker_nearby,",");
      String pos_str = pos2;
      int loc = tmp.Find(pos2);
      if(loc==-1)
	 return result;
      StringArray covs;
      covs.AddTokens(covReader.marker_cov,",");
      result = covs[loc].AsDouble()*SampleSize[study];
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
	 int loc = tmp.Find(pos);
	 if(loc==-1)
	    return result;
	 tmp.Clear();
	 tmp.AddTokens(covReader.markerNearbyCov[marker_idx],",");
	 result = tmp[loc].AsDouble()*SampleSize(study);
	 return result;
      }
      if(!covReader.ReadRecord(chr2,pos2.AsInteger()))
	 return result;
      StringArray tmp;
      tmp.AddTokens(covReader.marker_nearby,",");
      String pos_str = pos1;
      int loc = tmp.Find(pos1);
      if(loc==-1)
	 return result;
      StringArray covs;
      covs.AddTokens(covReader.marker_cov,",");
      result = covs[loc].AsDouble()*SampleSize[study];
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
   result.Dimension(commonVar_study[study].Length(),commonVar_study[study].Length(),0.0);
   for(int i=0;i<result.rows;i++)
   {
      result[i][i] = commonVar_marker_cov[study][i][0]*SampleSize[study];
      for(int j=i+1;j<result.cols;j++)
      {
	 if(common_chr[commonVar_study[study][j]] != common_chr[commonVar_study[study][i]])
	    continue;
	 if(common_pos[commonVar_study[study][i]]<common_pos[commonVar_study[study][j]])
	 {
	    int loc = commonVar_markers_in_window[study][i].Find(common_chr[commonVar_study[study][j]]);
	    if(loc==-1)
	       continue;
	    result[i][j] = result[j][i] = commonVar_marker_cov[study][i][loc]*SampleSize[study];
	 }
	 else
	 {
	    int loc = commonVar_markers_in_window[study][j].Find(common_chr[commonVar_study[study][i]]);
	    if(loc==-1)
	       continue;
	    result[i][j] = result[j][i] = commonVar_marker_cov[study][j][loc]*SampleSize[study];
	 }
      }
   }
   SVD svd;
   svd.InvertInPlace(result);
}

