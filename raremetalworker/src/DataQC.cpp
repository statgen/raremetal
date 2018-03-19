/**
 * Perform basic data quality checks. These include checking for consistency between the DAT and PED files,
 */

#include "DataQC.h"
#include "StringBasics.h"
#include "Error.h"
#include "PreMeta.h"

SanityCheck::SanityCheck()
{}

void SanityCheck::Check(Pedigree &ped, FILE *log)
{
    StringArray traits;
    if (!FastFit::traitName.IsEmpty())
    {
        FILE *trFile = fopen(FastFit::traitName, "r");
        if (trFile != nullptr)
        {
            String buffer;
            while (!feof(trFile))
            {
                buffer.ReadLine(trFile);
                traits.AddTokens(buffer, SEPARATORS);
            }
            fclose(trFile);
        } else
        {
            traits.AddTokens(FastFit::traitName, "/");
        }

        //check if all trait names are included in PED/DAT file.
        for (int i = 0; i < traits.Length(); i++)
        {
            if (ped.traitNames.Find(traits[i]) == -1)
            {
                traits.Delete(i);
                printf("Trait name %s is not included in DAT file, therefore excluded from analysis.\n",
                       traits[i].c_str());
                fprintf(log, "Trait %s is not included in DAT file, therefore excluded from analysis.\n",
                        traits[i].c_str());
            }
            if (traits.Length() == 0)
            {
                fprintf(log,
                        "ERROR! There is no trait to be analyzed. Please make sure your DAT file has the trait names you entered.\n");
                error("There is no trait to be analyzed. Please make sure your DAT file has the trait names you entered.\n");
            }
        }
    }


    if (ped.markerCount > 0)
    {
        PreMeta::genoFromPed = true;
    }

    if (FastFit::makeResiduals)
    {
        FastFit::useCovariates = true;
    }

    //DAT file is checked in PreMeta.cpp

    //STEP 2: summarizing VCF file.


    int vcf_marker_num = 0;
    if (PreMeta::vcfInput != "")
    {
        savvy::indexed_reader reader(PreMeta::vcfInput.c_str(), {""}, savvy::fmt::allele);

        if (!reader.good())
        {
            fprintf(log,
                    "FATAL ERROR! Either VCF does not exist or VCF file has to be tabix indexed. Bgzip your vcf file, then use \"tabix -p vcf your.vcf.gz\" to index.\n");
            error("Either VCF does not exist or VCF file has to be tabix indexed. Bgzip your vcf file, then use \"tabix -p vcf your.vcf.gz\" to index.\n");
        }


        for (const std::string &chr : reader.chromosomes())
        {
            chromosomeVCF.Push(chr.c_str());
        }



        //check if VCF file is empty

        if (chromosomeVCF.Length())
        {
            reader.reset_region({chromosomeVCF[0].c_str()});
            savvy::variant<std::vector<float>> record;
            while (reader >> record)
            {
                vcf_marker_num++;
                break;
            }
        }

        if (vcf_marker_num > 0)
        {
            PreMeta::genoFromVCF = true;
        }
    }

    //Check if VCF has any autosomal markers
    if (chromosomeVCF.Length() == 1 && chromosomeVCF[0] == PreMeta::xLabel)
    {
        //FastFit::autoOnly=true;
    }

    //STEP 3:Decide where to get marker genotype
    if (vcf_marker_num == 0 && PreMeta::vcfInput != "" && PreMeta::genoFromPed)
    {
        printf("Your vcf file has no marker included. Marker information will be read from PED file.\n");
        fprintf(log, "Your vcf file has no marker included. Marker information will be read from PED file.\n");
    } else if (vcf_marker_num == 0 && !PreMeta::genoFromPed)
    {
        fprintf(log,
                "ERROR! Neither your PED file nor your VCF file has marker information included. Raremetalworker needs marker genotypes to proceed.\n");
        error("Neither your PED file nor your VCF file has marker information included. Raremetalworker needs marker genotypes to proceed.\n");
    }

    if (PreMeta::genoFromPed && PreMeta::genoFromVCF)
    {
        PreMeta::genoFromPed = false;
        printf("  WARNING: Both of your PED and VCF files have marker genotypes. Markers will be read from VCF file. If you need markers included in PED file to be analyzed also, please merge the markers into either one PED file or one VCF file.\n");
        fprintf(log,
                "  WARNING: Both of your PED and VCF files have marker genotypes. Markers will be read from VCF file. If you need markers included in PED file to be analyzed also, please merge the markers into either one PED file or one VCF file.\n");
    }

    if (PreMeta::dosage && !PreMeta::genoFromVCF)
    {
        fprintf(log,
                "ERROR! Dosage has to be read from vcf file. Please make sure your dosage data is saved in a VCF input file.\n");
        error("Dosage has to be read from vcf file. Please make sure your dosage data is saved in a VCF input file.\n");
    }

    //STEP 2: check if there are genotypes included in input files at all.
    //printf("Deciding whether to treat them as related or unrelated samples\n");
    if (FastTransform::readInEmp != "")
    {
        FastTransform::empKin = true;
    }

    if (!FastTransform::empKin && !FastTransform::pedKin && FastTransform::readInEmp == "")
    {
        // If the user has not specified any of the arguments related to kinship calculation, remove attempts to account for it
        FastFit::unrelated = true;
    }
    printf("To account for relatedness in the sample,\n");
    fprintf(log, "To account for relatedness in the sample,\n");

    if (FastFit::unrelated)
    {
        printf(" NOTHING WILL BE DONE\n\n");
        fprintf(log, " NOTHING WILL BE DONE.\n\n");
    } else if (FastTransform::readInEmp != "")
    {
        printf(" EMPIRICAL KINSHIP WILL BE RETRIEVED FROM FILE.\n\n");
        fprintf(log, " EMPIRICAL KINSHIP WILL BE RETRIEVED FROM FILE.\n\n");
    } else if (FastTransform::empKin)
    {
        printf(" EMPIRICAL KINSHIP WILL BE ESTIMATED.\n\n");
        fprintf(log, " EMPIRICAL KINSHIP WILL BE ESTIMATED.\n\n");
    } else if (FastTransform::pedKin)
    {
        printf(" PEDIGREE INFORMATION WILL BE USED.\n\n");
        fprintf(log, " PEDIGREE INFORMATION WILL BE USED.\n\n");
    }

    printf("Sanity checking input files ...\n");
    fprintf(log, "Sanity checking input files ...\n");

    if (FastTransform::mergedVCFID)
    {
        printf("  WARNING: --mergedVCFID has been activated; RAREMETALWORKER assumes sample IDs in VCF file are formated FAMID_PID.\n");
        fprintf(log,
                "  WARNING: --mergedVCFID has been activated; RAREMETALWORKER assumes sample IDs in VCF file are formated FAMID_PID.\n");
    }

    //STEP 3: Check if there are overlapping sample IDs from PED and VCF file
    if (PreMeta::genoFromVCF)
    {
        bool status = false;
        StringIntHash pedSampleID;
        for (int s = 0; s < ped.count; s++)
        {
            if (FastTransform::mergedVCFID)
            {
                pedSampleID.SetInteger(ped[s].famid + "_" + ped[s].pid, s);
            } else
            {
                pedSampleID.SetInteger(ped[s].pid, s);
            }
        }

        savvy::reader reader(PreMeta::vcfInput.c_str(), savvy::fmt::allele);
        int numSamples = reader.samples().size();
        for (int s = 0; s < numSamples; s++)
        {
            const char *sample = reader.samples()[s].c_str();
            if (pedSampleID.Integer(sample) >= 0)
            {
                status = true;
                break;
            }
        }

        if (!status)
        {
            fprintf(log,
                    "ERROR! There is no overlapping sample ID between PED and VCF file. Please double check to make sure that PED file and VCF file have some or all sample IDs matched.\n");
            error("There is no overlapping sample ID between PED and VCF file. Please double check to make sure that PED file and VCF file have some or all sample IDs matched.\n");
        }
    }

    //STEP 4:if X chromosome is included then check quality of genotype for chromosome X
    //printf("Checking the quality of genotype from chr X.\n");
    //bool Xstatus=false;
    //if(AutoFit::fitX || PreMeta::xLabel!="X")
    //	Xstatus=true;

    if (AutoFit::fitX && FastTransform::readInEmp != "")
    {
        if (FastTransform::readInEmpX == "")
        {
            //error("--vcX and --kinFile options are used, but --kinxFile option is not specified. Please enter your kinshipX file following --kinxFile, or remove --vcX to make RAREMETAL run. Another solution is to estimate kinship matrices using --kinGeno --kinSave --vcX options.");
        }
    }

    if (FastTransform::readInEmpX != "" && !AutoFit::fitX)
    {
        error("  WARNING: --kinxFile has to be used with --vcX option. --kinxFile option will be ignored.\n");
        //FastTransform::readInEmpX = "";
    }

    if (AutoFit::fitX && FastTransform::readInEmpX == "" && FastTransform::readInEmp != "")
    {
        String tmp = FastTransform::readInEmp;
        int loc = tmp.Find("Kinship");
        FastTransform::readInEmpX = tmp.SubStr(0, loc) + "KinshipX";
        printf("  WARNING: --kinxFile is not specified. Will read kinshipX file from %s. Or, please enter kinshipX file in --kinxFile.\n",
               FastTransform::readInEmpX.c_str());
        fprintf(log,
                "  WARNING: --kinxFile is not specified. Will read kinshipX file from %s. Or, please enter kinshipX file in --kinxFile.\n",
                FastTransform::readInEmpX.c_str());
    }

    if (FastFit::separateX && !AutoFit::fitX)
    {
        error("WARNING: --separateX option has to be used with --vcX. --separateX option will be deactiaved.\n");
        //FastFit::separateX=false;
    }

    if (FastFit::separateX)
    {
        AutoFit::fitX = false;
    }

    //part 1: check if label for X chromosome is correct
    if (PreMeta::genoFromVCF)
    {
        int xpos = -1;
        for (int i = 0; i < chromosomeVCF.Length(); i++)
        {
            if (chromosomeVCF[i].FindChar('X') != -1 || chromosomeVCF[i].FindChar('x') != -1)
            {
                xpos = i;
                //		Xstatus=true;
                break;
            }
        }

        if (chromosomeVCF.Find(PreMeta::xLabel) == -1 && xpos != -1)
        {
            fprintf(log,
                    "ERROR! Please double check label for X chromosome. Your VCF has X chromosome labeled as %s, but your --xLabel option has entry %s.\n",
                    chromosomeVCF[xpos].c_str(), PreMeta::xLabel.c_str());
            error("Please double check label for X chromosome. Your VCF has X chromosome labeled as %s, but your --xLabel option has entry %s.\n",
                  chromosomeVCF[xpos].c_str(), PreMeta::xLabel.c_str());
        }
    }

    if (PreMeta::genoFromPed)
    {
        String xname = "";
        for (int i = 0; i < ped.markerCount; i++)
        {
            if (ped.markerNames[i].FindChar('X') != -1 || ped.markerNames[i].FindChar('x') != -1)
            {
                //		Xstatus=true;
                StringArray name;
                name.AddTokens(ped.markerNames[i], ":");
                xname = name[0];
                break;
            }
        }

        if (xname != PreMeta::xLabel && xname != "")
        {
            fprintf(log,
                    "ERROR! Please double check label for X chromosome. Your VCF has X chromosome labeled as %s, but your --xLabel option has entry %s.\n",
                    xname.c_str(), PreMeta::xLabel.c_str());
            error("Please double check label for X chromosome. Your VCF has X chromosome labeled as %s, but your --xLabel option has entry %s.\n",
                  xname.c_str(), PreMeta::xLabel.c_str());
        }
    }

    //part 2: if --vcX option specified or --xLabel specified, then check to make sure there are markers on chromosome X.

    if (AutoFit::fitX || PreMeta::xLabel != "X")
    {
        if (PreMeta::genoFromVCF)
        {
            if (chromosomeVCF.Find(PreMeta::xLabel) == -1)
            {
                fprintf(log,
                        "  WARNING: --vcX was requested, but no marker from chromosome X was labeled \"%s\" in your VCF file.\n",
                        PreMeta::xLabel.c_str());
                printf("  WARNING: --vcX was requested, but no marker from chromosome X was labeled \"%s\" in your VCF file.\n",
                       PreMeta::xLabel.c_str());
            }
        } else
        {
            StringArray markernames;
            bool status = false;
            for (int m = 0; m < ped.markerCount; m++)
            {
                markernames.Clear();
                markernames.AddTokens(ped.markerNames[m], ":");
                if (markernames[0] == PreMeta::xLabel)
                {
                    status = true;
                    break;
                }
            }
            if (!status)
            {
                fprintf(log,
                        "  WARNING: --vcX was requested, but no marker from chromosome X was labeled \"%s\" in your input file.\n",
                        PreMeta::xLabel.c_str());
                printf("  WARNING: --vcX was requested, but no marker from chromosome X was labeled \"%s\" in your input file.\n",
                       PreMeta::xLabel.c_str());
            }
        }
    }

    //part 2: if male X chromosome has genotype not in 0/0 or 1/1 or ./. or ./0 or ./1 format then warning.
    //check heterozygosity among males. If over 10% then ERROR. Otherwise skip het ones.

    /*This is for screening out bad nonPAR variants in male
      if(Xstatus)
      {
      int bad = 0,total=0;
      if(PreMeta::genoFromVCF)
      {
    //get all male sample from PED file first.
    StringIntHash maleID;
    for(int i=0;i<ped.count;i++)
    {
    if(ped[i].sex==PreMeta::maleLabel)
    maleID.SetInteger(ped[i].pid,i);
    }
    reader.open(PreMeta::vcfInput,header);
    int numSamples = header.getNumSamples();
    while(reader.readRecord(record))
    {
    String chr = record.getChromStr();
    int pos = record.get1BasedPosition();
    if(chr!=PreMeta::xLabel)
    continue;
    if(pos<PreMeta::Xstart || pos>PreMeta::Xend)
    continue;
    total++;
    bool badsnp=false;
    for(int s=0;s<numSamples;s++)
    {
    const char * sample = header.getSampleName(s);
    int numGTs = record.getNumGTs(s);
    bool missing =false;
    if(maleID.Integer(sample)>=0)
    {
    int sum = 0;
    for(int j = 0; j < numGTs; j++)
    {
    int a = record.getGT(s,j);
    if(a==VcfGenotypeSample::MISSING_GT)
    {
    sum=0;
    missing =true;
    break;
    }
    if(a>1)
    {
    badsnp=true;
    bad++;
    printf("WARNING: variant %s:%d has allele %d and will skipped from analysis.\n",chr.c_str(),pos);
    fprintf(log,"WARNING: variant %s:%d has allele %d and will skipped from analysis.\n",chr.c_str(),pos);
    break;
    }
    sum += a;
    }
    //if at least one genotype is 1/1 then bad++
    //if((!missing && numGTs==2 && sum%2!=0) || badsnp)
    if(badsnp)
    {
    String SNPname;
    SNPname =  chr + ":" +pos;
    skippedSNPs.SetInteger(SNPname,bad);
    break;
    }
    }
    }
    }
    reader.close();
    }

    if(PreMeta::genoFromPed)
    {
    StringArray markername;
    String chr;
    int pos;
    for(int m=0;m<ped.markerCount;m++)
    {
       markername.AddTokens(ped.markerNames[m],":");
       chr = markername[0];
       pos = atoi(markername[1].c_str());
       if(chr!=PreMeta::xLabel)
      continue;
       if(pos<PreMeta::Xstart || pos>PreMeta::Xend)
      continue;
       total++;
       for(int i=0;i<ped.count;i++)
       {
      if(ped[i].sex==PreMeta::maleLabel)
      {
         if(!ped[i].isGenotyped(m))
            continue;
         if(ped[i].markers[m].isHeterozygous())
         {
            bad++;
            break;
         }
      }
       }
    }
 }
 if(total>0)
 {
    if(bad>total*0.2)
    {
       fprintf(log,"FATAL ERROR: There are over 20% variants on non-pseudo-autosomal region of chromosome X with heterozygous male. Genotype for male samples in this region should be ./0, ./1, ./., 0, 1, 1/1 or 0/0. Analsysis will abort.\n");
       error("There are over 20% variants on non-pseudo-autosomal region of chromosome X with heterozygous male. Genotype for male samples in this region should be ./0, ./1, ./., 0, 1, 1/1 or 0/0. Analsysis will abort.\n");
    }
 }

 if(bad>0)
 {
    fprintf(log,"\tWARNING: There are %d variants out of %d variants on non-pseudo-autosomal region of chromosome X with heterozygous male genotype. Genotype for male samples in this region should be ./0, ./1, ./., 0, 1, 1/1 or 0/0. These variants will be skipped.\n",bad,total);
    printf("\tWARNING: There are %d variants out of %d variants on non-pseudo-autosomal region of chromosome X with heterozygous male genotype. Genotype for male samples in this region should be ./0, ./1, ./., 0, 1, 1/1 or 0/0. These variants will be skipped.\n",bad,total);
 }
 }
 */
    printf("completed.\n\n");
    fprintf(log, "completed.\n\n");
}
