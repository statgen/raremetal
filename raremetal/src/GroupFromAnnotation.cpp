#include "GroupFromAnnotation.h"
#include "Error.h"
#include "MathVector.h"
#include "IntArray.h"
#include "InputFile.h"
#include "QuickIndex.h"

String GroupFromAnnotation::groupFile = "";
String GroupFromAnnotation::vcfInput = "";
String GroupFromAnnotation::function = "";
String GroupFromAnnotation::mapFile = "../data/refFlat_hg19.txt";
bool GroupFromAnnotation::labelHits = false;

GroupFromAnnotation::GroupFromAnnotation()
{
	geneCount = 0;
}

GroupFromAnnotation::~GroupFromAnnotation()
{
   //   if(SNPlist) delete [] SNPlist;
   //   if(SNPNoAllele) delete [] SNPNoAllele;
}

void GroupFromAnnotation::GetGeneMap(String path)
{
   IFILE genemap;
   genemap =  ifopen(mapFile,"r");
   if(genemap==NULL)
   {
      if(mapFile=="../data/refFlat_hg19.txt")
      {
	 mapFile += ".gz";

	 genemap = ifopen(mapFile,"r");
	 if(genemap==NULL)
	 {
	    int loc = path.Find("bin");
	    if(loc!=-1)
	    {
	       mapFile = path.Left(loc-1);
	       mapFile += "/data/refFlat_hg19.txt";
	    }
	    else
	    {
	       mapFile += "../data/refFlat_hg19.txt";
	    }
	    genemap = ifopen(mapFile,"r");
	 }
	 if(genemap==NULL)
	 {
	    mapFile += ".gz";
	    genemap = ifopen(mapFile,"r");
	 }
	 if(genemap==NULL)
	    error("Cannot open gene mapping file %s.\n",mapFile.c_str());

      }
      else
	 error("Cannot open gene mapping file %s.\n",mapFile.c_str());
   }
   StringIntHash GeneLocHash;
   StringArray strand;
   int gene_idx =0;

   while(!ifeof(genemap))
   {
      String buffer;
      buffer.ReadLine(genemap);
      StringArray record;
      record.AddTokens(buffer,"\t");
      int loc = GeneLocHash.Integer(record[0]);
      if(loc==-1)
      {
	 GeneLocHash.SetInteger(record[0],gene_idx);
	 //save chr, start and end positions
	 StringArray gene_chr;
	 if(record[2][2]=='r' || record[2][2]=='R')
	    record[2] = record[2].SubStr(3);
	 gene_chr.AddTokens(record[2],"_,;.");
	 if(gene_chr[0].Find("Un")!=-1)
	    continue;
	 /*
	    if(ChrLocHash.Integer(gene_chr[0])==-1)
	    {
	    chr_count++;
	    unique_chr.Push(gene_chr[0]);
	    ChrLocHash.SetInteger(gene_chr[0],chr_count);
	    }
	  */
	 chr.Push(gene_chr[0]);
	 //printf("%d\t%s\t%s\n",idx,record[0].c_str(),gene_chr[0].c_str());
	 start_pos.Push(record[4].AsInteger());
	 end_pos.Push(record[5].AsInteger());
	 strand.Push(record[3]);
	 genename.Push(record[0]);
	 gene_idx++;
      }
      else
      {
	 //get the current chr
	 StringArray gene_chr;
	 if(record[2][2]=='r' || record[2][2]=='R')
	    record[2] = record[2].SubStr(3);
	 gene_chr.AddTokens(record[2],"_,;.");
	 if(gene_chr[0].Find("Un")!=-1)
	    continue;
	 //check if strand and chr are consistent with previous record
	 if(chr[loc]!=gene_chr[0]) 
	    //if(strand[loc]!=record[3] || chr[loc]!=gene_chr[0]) 
	    //    printf("Gene %s in %s has multiple records in different chromosome or strand.\n",record[0].c_str(),mapFile.c_str());
	    continue;
	 //update start and end position
	 if(record[4].AsInteger()<start_pos[loc])
	    start_pos[loc] = record[4].AsInteger();
	 if(record[5].AsInteger()>end_pos[loc])
	    end_pos[loc] = record[5].AsInteger();
      }
   }
   ifclose(genemap);
   //ifclose(genemap);
   chr_idx.Index(chr);
   String chr_=chr[chr_idx[0]];
   for(int i=1;i<chr.Length();i++)
   {
      if(chr[chr_idx[i]]!=chr_)
      {
	 ChrStartHash.SetInteger(chr[chr_idx[i]],i);
	 ChrEndHash.SetInteger(chr_,i-1);
	 chr_ = chr[chr_idx[i]];
      }
   }
}

String GroupFromAnnotation::AnnotateSingleVar(String chr, int pos)
{
   String gene_belonged = "";
   int start = ChrStartHash.Integer(chr);
   int end = ChrEndHash.Integer(chr);
   for(int i=start;i<=end;i++)
   {
      int idx = chr_idx[i];
      if(pos>=start_pos[idx] && pos<=end_pos[idx])
      {
	 if(gene_belonged!="")
	 {
	    //  if(gene_belonged.Find('/')!=-1)
	    continue;
	    //gene_belonged+="/";
	    //gene_belonged += group.genename[idx];
	 }
	 else
	    gene_belonged += genename[idx];
      }
   }
   return gene_belonged;
}

void GroupFromAnnotation::GetGroupFromFile(FILE * log)
{
   //Fill in annoGroups.
   StringArray tmp;
   FILE * file = fopen(groupFile,"r");
   if(file==NULL)
   {
      printf("ERROR! Cannot open group file %s.\n",groupFile.c_str());
      error("ERROR! Cannot open group file %s.\n",groupFile.c_str());
   }
   String buffer;
   int line = 0;
   while (!feof(file))
   {
      buffer.ReadLine(file);
      tmp.Clear();
      tmp.AddTokens(buffer, SEPARATORS);
      if(tmp.Length()==0)
	 continue;
      annoGroups.Push(tmp[0]);
      chrom.Push(tmp[1]);
      line++;
   }
   fclose(file);

   //Fill in SNPlist.
   SNPlist = new StringArray [line];
   SNPNoAllele = new StringArray [line];
   FILE * samefile = fopen(groupFile,"r");
   line = 0;
   Vector pos;
   while (!feof(samefile))
   {
      buffer.ReadLine(samefile);
      tmp.Clear();
      pos.Clear();
      tmp.AddTokens(buffer, "\t ");
      SNPlist[line].Dimension(0);
      SNPNoAllele[line].Dimension(0);
      for(int i=1;i<tmp.Length();i++)
      {
	 SNPlist[line].Push(tmp[i]);
	 StringArray sub;
	 sub.Clear();
	 sub.AddTokens(tmp[i],":_/");
	 if(sub.Length()!=4)
	 {
	    printf("Warning: group %s has a variant %s that has invalid format. The correct format should be chr:pos:allele1:allele2.\n",tmp[0].c_str(),tmp[i].c_str());
	    fprintf(log,"Warning: group %s has a variant %s that has invalid format. The correct format should be chr:pos:allele1:allele2.\n",tmp[0].c_str(),tmp[i].c_str());
	    continue;
	 }
	 pos.Push(sub[1].AsInteger());
	 SNPNoAllele[line].Push(sub[0] + ":" + sub[1]);
      }
      //sort SNPlist[line] and SNPNoAllele[line]
      if(SNPlist[line].Length()>1)
      {
	 Vector sorted_pos,order;
	 sorted_pos.Copy(pos);
	 sorted_pos.Sort();
	 order.Dimension(pos.Length());
	 for(int i=0;i<sorted_pos.Length();i++)
	 {
	    for(int j=0;j<pos.Length();j++)
	    {
	       if(sorted_pos[i]==pos[j])
	       {
		  order[i]=j; 
		  break;
	       }
	    }
	 }

	 StringArray cp_SNPlist,cp_SNPNoAllele;
	 cp_SNPlist.Dimension(SNPlist[line].Length());
	 cp_SNPNoAllele.Dimension(SNPNoAllele[line].Length());
	 for(int l=0;l<SNPlist[line].Length();l++)
	 {
	    cp_SNPlist[l] = SNPlist[line][l];
	    cp_SNPNoAllele[l] = SNPNoAllele[line][l];
	 }
	 for(int i=0;i<order.Length();i++)
	 {
	    SNPlist[line][i] = cp_SNPlist[order[i]];
	    //printf("%s\t",SNPlist[line][i].c_str());
	    SNPNoAllele[line][i] = cp_SNPNoAllele[order[i]] ;
	 }
	 //printf("\n");
      }
      line++;
   }
   fclose(samefile);
}

/* Original function
void GroupFromAnnotation::GetGroupFromVCF()
{  
   printf("Parsing annotations from annotated VCF file ...\n");
   StringArray func;
   func.AddTokens(function,"/");

   FILE * inFile;
   inFile = fopen(vcfInput,"r");
   String buffer;
   StringArray tmp;
   StringIntHash groupHash;
   int geneCount=0;

   while (!feof(inFile))
   {
      buffer.ReadLine(inFile);
      tmp.Clear();
      if(function != "")
      {
	 tmp.AddTokens(buffer, "\t;/= ");
	 for(int i=0;i<tmp.Length();i++)
	 {
	    int result = -1;
	    for(int f =0;f<func.Length();f++)
	    {
	       String tmp_i,tmp_i_upper,func_f_upper;
	       tmp_i = tmp[i].SubStr(0,func[f].Length());
	       tmp_i_upper = tmp_i.ToUpper();
	       func_f_upper = func[f].ToUpper();
	       result = tmp_i.SlowFind(func_f_upper);
	       if(result==0)
		  break;
	    }
	    if(result == 0)
	    {
	       StringArray sub;
	       sub.AddTokens(tmp[i],":(,");
	       String genename = sub[1];
	       int gene = groupHash.Integer(genename);
	       if(gene==-1)
	       {
		  groupHash.SetInteger(genename,geneCount);
		  geneCount++;
		  annoGroups.Push(genename);
		  chrom.Push(tmp[0]);
	       }
	    }
	 }
      }
      else
      {
	 tmp.AddTokens(buffer, "\t;/ ");
	 for(int i=0;i<tmp.Length();i++)
	 {
	    StringArray tmp_i;
	    tmp_i.AddTokens(tmp[i],"=:(,");
	    if(tmp_i[0]=="ANNO")
	    {
	       String genename = tmp_i[2];
	       int gene = groupHash.Integer(genename);
	       if(gene==-1)
	       {
		  groupHash.SetInteger(genename,geneCount);
		  geneCount++;
		  annoGroups.Push(genename);
		  chrom.Push(tmp[0]);
	       }
	    }
	 }
      }
   }
   fclose(inFile);

   inFile = fopen(vcfInput,"r");
   SNPlist = new StringArray [geneCount];
   SNPNoAllele = new StringArray [geneCount];
   Vector * pos;
   pos = new Vector [geneCount];

   while (!feof(inFile))
   {
      buffer.ReadLine(inFile);
      tmp.Clear();
      if(function != "")
      {
	 tmp.AddTokens(buffer, "\t;/= ");
	 int oldgene = -1;
	 for(int i=0;i<tmp.Length();i++)
	 {
	    int result = -1;
	    for(int f =0;f<func.Length();f++)
	    {
	       String tmp_i;
	       tmp_i = tmp[i].SubStr(0,func[f].Length());
	       result = tmp_i.SlowFind(func[f]);
	       if(result==0)
		  break;
	    }

	    if(result==0)
	    {
	       StringArray sub;
	       sub.AddTokens(tmp[i],":(,");
	       String genename = sub[1];
	       int gene = groupHash.Integer(genename);
	       if(gene==oldgene)
		  continue;
	       oldgene =gene;

	       String markername;
	       if(tmp[0].Find("chr")!=-1)
		  markername = tmp[0].SubStr(3);
	       else
		  markername = tmp[0];

	       markername += ":";
	       markername += tmp[1];
	       markername += ":";
	       markername += tmp[3];
	       markername += ":";
	       markername += tmp[4];
	       SNPlist[gene].Push(markername);
	       SNPNoAllele[gene].Push(tmp[0]+":"+tmp[1]);
	       pos[gene].Push(tmp[1].AsInteger());
	    }
	 }
      }
      else
      {
	 tmp.AddTokens(buffer, "\t;/ ");
	 int oldgene = -1;
	 for(int i=0;i<tmp.Length();i++)
	 {
	    StringArray tmp_i;
	    tmp_i.AddTokens(tmp[i],"=:(,");
	    if(tmp_i[0]=="ANNO")
	    {
	       String genename = tmp_i[2];
	       int gene = groupHash.Integer(genename);
	       if(gene==oldgene)
		  continue;
	       oldgene =gene;

	       String markername;
	       if(tmp[0].Find("chr")!=-1)
		  markername = tmp[0].SubStr(3);
	       else
		  markername = tmp[0];

	       markername += ":";
	       markername += tmp[1];
	       markername += ":";
	       markername += tmp[3];
	       markername += ":";
	       markername += tmp[4];
	       SNPlist[gene].Push(markername);
	       SNPNoAllele[gene].Push(tmp[0]+":"+tmp[1]);
	       pos[gene].Push(tmp[1].AsInteger());
	    }
	 }
      }
   }
   fclose(inFile);

   //sort SNPlist and SNPNoAllele
   for(int line=0;line<geneCount;line++)
   {
      if(SNPlist[line].Length()>1)
      {
	 Vector sorted_pos,order;
	 sorted_pos.Copy(pos[line]);
	 sorted_pos.Sort();
	 order.Dimension(pos[line].Length());
	 for(int i=0;i<sorted_pos.Length();i++)
	 {
	    for(int j=0;j<pos[line].Length();j++)
	    {
	       if(sorted_pos[i]==pos[line][j])
	       {
		  order[i]=j;
		  break;
	       }
	    }
	 }

	 StringArray cp_SNPlist,cp_SNPNoAllele;
	 cp_SNPlist.Dimension(SNPlist[line].Length());
	 cp_SNPNoAllele.Dimension(SNPNoAllele[line].Length());
	 for(int l=0;l<SNPlist[line].Length();l++)
	 {
	    cp_SNPlist[l] = SNPlist[line][l];
	    cp_SNPNoAllele[l] = SNPNoAllele[line][l];
	 }
	 for(int i=0;i<order.Length();i++)
	 {
	    SNPlist[line][i] = cp_SNPlist[order[i]];
	    SNPNoAllele[line][i] = cp_SNPNoAllele[order[i]] ;
	 }
      }
   }
   printf("done!\n");

  //  This is for outpupt groupfile from annotated vcf file 
      FILE * out;
      String filename = "test.groupfile";
      out = fopen(filename,"wt");
      for(int i=0;i<geneCount;i++)
      {
      fprintf(out,"%s\t",annoGroups[i].c_str());
      for(int m=0;m<SNPlist[i].Length();m++)
      fprintf(out,"%s\t",SNPlist[i][m].c_str());
      fprintf(out,"\n");
      }
      fclose(out);
    
}
*/


// partition 8th column only
void GroupFromAnnotation::GetGroupFromVCF()
{  
   printf("Parsing annotations from annotated VCF file ...\n");
   StringArray func;
   func.AddTokens(function,"/");
   
   vcfInitialize(); // set size of the tables

   FILE * inFile;
   inFile = fopen(vcfInput,"r");
   StringIntHash groupHash;
   int geneCount=0;

// add all genes to group hash first
	while (!feof(inFile))
	{
		String buffer;
		buffer.ReadLine(inFile);
		if ( buffer[0] == '#' )
			continue;
		addLineFromVcf( buffer );
	}
	fclose(inFile);
	
// sort SNPlist and SNPNoAllele
	for( int g=0; g<geneCount; g++ ) {
		if ( SNPlist[g].Length()>1 ) {
			Vector order;
			setOrderFromSortedPositions( g, order);
			
			StringArray cp_SNPlist,cp_SNPNoAllele;
			cp_SNPlist.Dimension(SNPlist[g].Length());
			cp_SNPNoAllele.Dimension(SNPNoAllele[g].Length());
			for(int l=0;l<SNPlist[g].Length();l++) {
				cp_SNPlist[l] = SNPlist[g][l];
				cp_SNPNoAllele[l] = SNPNoAllele[g][l];
			}
			for(int i=0;i<order.Length();i++) {
				SNPlist[g][i] = cp_SNPlist[order[i]];
				SNPNoAllele[g][i] = cp_SNPNoAllele[order[i]] ;
			}
		}
	}

// print test.groupfile
	printf("done!\n");
	String grp_filename = "test.groupfile";
	printGroupFile( grp_filename );
}	


void GroupFromAnnotation::addLineFromVcf( String & buffer )
{
// sample:
// ANNO=Nonsynonymous:ASB16;
// ANNOFULL=ASB16/NM_080863:+:Nonsynonymous(CCC/Pro/P->ACC/Thr/T:Base1310/1362:Codon437/454:Exon5/5):Exon
//	|C17orf65/NM_178542:-:Intron

	StringArray vfield;
	vfield.AddTokens(buffer, "\t");
	if ( vfield.Length() < 8 )
		error("Annotation vcf only has %d columns!\n", vfield.Length());
	StringArray info_semicolon;
	info_semicolon.AddTokens( vfield[7],";" );
	
// find ANNOFULL first
	int annofull_index = -1;
	for( int i=0; i<info_semicolon.Length(); i++ ) {
		String iheader = info_semicolon[i].SubStr(0,8);
		if (iheader == "ANNOFULL") {
			annofull_index = i;
			break;
		}
	}
	if (annofull_index == -1) {
		printf("warning: no ANNOFULL field at chr%s:%s. Variant won't included in groups!\n", info_semicolon[0].c_str(), info_semicolon[1].c_str());
		return;
	}

// remove ANNOFULL=
	String anno_full_str = info_semicolon[annofull_index].SubStr(9);

// check each alternative field
	StringArray alts;
	alts.AddTokens( anno_full_str, "|" );
	for( int a=0; a<alts.Length(); a++ ) {
		StringArray sub;
		sub.AddTokens( alts[a], ":/=");
		if (func_upper.Length() != 0) { // match before add
			for(int f =0;f<func_upper.Length();f++) {
				bool pattern_match = checkPatternMatch( sub, func_upper[f] );
				if ( pattern_match ) {
					addGeneFromVcf( vfield, sub[0] );
					break;
				}
			}
		}
		else { // no pattern to match: check if intergenic first
			String upper_name = sub[0].ToUpper();
			if ( !upper_name.SlowFind( "INTERGENIC" ) )
				addGeneFromVcf( vfield, sub[0] );
		}
	}	
}

bool GroupFromAnnotation::checkPatternMatch( StringArray & sub, String & func )
{
	int result = -1;
	for( int i=0; i<sub.Length(); i++ ) {
		if ( sub[i].Length() < func.Length() )
			continue;
		String str = sub[i];
		String upper_sub = str.ToUpper();
		result = upper_sub.SlowFind( func );
		if ( result == 0 )
			break;
	}
	if ( result == 0 )
		return 1;
	else
		return 0;
}
		
void GroupFromAnnotation::addGeneFromVcf( StringArray & vfield, String & genename )
{
	int gene = groupHash.Integer( genename );
	if ( gene == -1 )
		error("Gene %s not recorded..something is wrong!\n", genename.c_str());
	String markername;
	if (vfield[0].Find("chr") != -1)
		markername = vfield[0].SubStr(3);
	else
		markername = vfield[0];
	markername += ":";
	markername += vfield[1];
	markername += ":";
	markername += vfield[3];
	markername += ":";
	markername += vfield[4];
	if ( SNPlist[gene].FastFind(markername) == -1 ) { // do not add duplicate
		SNPlist[gene].Push(markername);
		SNPNoAllele[gene].Push(vfield[0]+":"+vfield[1]);
		pos[gene].Push(vfield[1].AsInteger());
	}
}		

void GroupFromAnnotation::addGeneToGroupHash( String & gene)
{
	int gidx = groupHash.Integer( gene);
	if (gidx != -1)
		return;
	groupHash.SetInteger(gene, geneCount);
	annoGroups.Push(gene);
	geneCount++;
}

void GroupFromAnnotation::setOrderFromSortedPositions( int pos_idx, Vector & order )
{
	Vector sorted_pos;
	sorted_pos.Copy( pos[pos_idx] );
	sorted_pos.Sort();

	order.Dimension( sorted_pos.Length() );
	for( int i=0; i<sorted_pos.Length(); i++ ) {
		for( int j=0; j<sorted_pos.Length(); j++ ) {
			if ( sorted_pos[i] == pos[pos_idx][j] ) {
				order[i] = j;
				break;
			}
		}
	}
}


void GroupFromAnnotation::printGroupFile( String & filename )
{
	FILE * out;
    out = fopen(filename,"wt");
    for(int i=0;i<geneCount;i++) {
		fprintf(out,"%s\t",annoGroups[i].c_str());
		for(int m=0;m<SNPlist[i].Length();m++)
		fprintf(out,"%s\t",SNPlist[i][m].c_str());
		fprintf(out,"\n");
    }
    fclose(out);
}	

void GroupFromAnnotation::vcfInitialize()
{
	// func_upper
	if ( function != "" ) {
		func_upper.AddTokens( function, "/" );
		for( int i=0; i<func_upper.Length(); i++ )
			func_upper[i] = func_upper[i].ToUpper();
	}

	FILE * inFile;
	inFile = fopen(vcfInput,"r");
	while (!feof(inFile)) {
		String buffer;
		buffer.ReadLine( inFile);
		if ( buffer[0] == '#' )
			continue;
		StringArray vfield;
		vfield.AddTokens(buffer, "\t");
		if ( vfield.Length() < 8 )
			error("Annotation vcf only has %d columns!\n", vfield.Length());
		StringArray info_semicolon;
		info_semicolon.AddTokens( vfield[7],";" );
		
		int annofull_index = -1;
		for( int i=0; i<info_semicolon.Length(); i++ ) {
			String iheader = info_semicolon[i].SubStr(0,8);
			if (iheader == "ANNOFULL") {
				annofull_index = i;
				break;
			}
		}
		if (annofull_index == -1)
			continue;
		String anno_full_str = info_semicolon[annofull_index].SubStr(9);
		StringArray alts;
		alts.AddTokens( anno_full_str, "|" );
		for( int a=0; a<alts.Length(); a++ ) {
			StringArray sub;
			sub.AddTokens( alts[a], ":/=");
			if (func_upper.Length() != 0) { // match before add
				for(int f =0;f<func_upper.Length();f++) {
					bool pattern_match = checkPatternMatch( sub, func_upper[f] );
					if ( pattern_match ) {
						chrom.Push( vfield[0] );
						addGeneToGroupHash( sub[0] );
						break;
					}
				}
			}
			else { // no pattern to match
				chrom.Push( vfield[0] );
				addGeneToGroupHash( sub[0] );		
			}
		}
	}

// vectors	
	SNPlist = new StringArray [geneCount];
	SNPNoAllele = new StringArray [geneCount];
	pos = new Vector [geneCount];
	
}

void GroupFromAnnotation::Run(String path,FILE * log)
{ 
   if(vcfInput !="")
      GetGroupFromVCF();
   if(groupFile!="")
      GetGroupFromFile(log);
   if(labelHits)
      GetGeneMap(path);
}


