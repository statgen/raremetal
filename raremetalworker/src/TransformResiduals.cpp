////////////////////////////////////////////////////////////////////// 
// TransformResiduals.cpp 
// (c) 2012-2013 Shuang Feng, Dajiang Liu, Goncalo Abecasis
// 
// This file is distributed as part of the RareMetalWorker source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile RareMetalWorker.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday November 28, 2012
// 

#include "TransformResiduals.h"
#include "MathSVD.h"
#include "OutputKin.h"
#include "PreMeta.h"
#include "KinshipEmp.h"
#include <Eigen/Eigenvalues>
#include <map>
#include <iterator>

bool FastTransform::empKin = false;
bool FastTransform::pedKin = false;
bool FastTransform::mergedVCFID = false;
String FastTransform::readInEmp = "";
String FastTransform::readInEmpX = "";

FastTransform::FastTransform()
        : pPheno(NULL)
{}

FastTransform::~FastTransform()
{
    if (pPheno)
    { delete[] pPheno; }
}

void FastTransform::GetFixedEffectsCount(Pedigree &ped, bool useCovariates)
{
    fixedEffects = 1 + (useCovariates ? ped.covariateCount : 0);
}

void FastTransform::SelectSamplesPED(Pedigree &ped, bool useCovariates)
{
    //printf("Selecting useful samples ...\n");
    for (int i = 0; i < ped.count; i++)
    {
        //check if this sample is phenotyped at least one trait
        bool phenotyped = false;
        for (int tr = 0; tr < ped.traitNames.Length(); tr++)
        {
            if (ped[i].isPhenotyped(tr))
            {
                phenotyped = true;
                break;
            }
        }
        if (phenotyped)
        {
            bool reject = true;
            for (int m = 0; m < ped.markerCount; m++)
            {
                if (ped[i].isGenotyped(m))
                {
                    reject = false;
                    break;
                }
            }
            if (reject)
            {
                continue;
            }
            //String id = ped[i].famid + "." + ped[i].pid;
            //samplePEDIDHash.SetInteger(id,i);
            samplePEDIDHash.SetInteger(ped[i].pid, i);
            genotypedSamplePED.Push(i);
        }
    }
}


// 07/08/16 update: read vcf once instead of multiple times
void FastTransform::SelectSamplesVCF(Pedigree &ped, bool useCovariates)
{
    // printf("Selecting useful samples ...\n");


    savvy::reader reader(PreMeta::vcfInput.c_str(), PreMeta::dosage ? savvy::fmt::dosage : savvy::fmt::genotype);
    int numSamples = reader.samples().size();
    totalN = numSamples;

    StringIntHash VCFID;
    for (int s = 0; s < numSamples; s++)
    {
        const char *sample = reader.samples()[s].c_str();
        VCFID.SetInteger(sample, s);
    }

    // now see if each sample is genotyped in at least one site
    std::map<int, bool> g1; // store the samples that haven't found a genotyped site
    for (int s = 0; s < numSamples; s++)
    {
        g1[s] = 0;
    }

    savvy::site_info record_anno;
    std::vector<float> record_data;
    while (reader.read(record_anno, record_data))
    {
        if (g1.empty())
        {
            break;
        }

        if (record_data.empty() && PreMeta::dosage)
        {
            error("cannot find dosage at \"%s\" in VCF!", PreMeta::dosageFlag.c_str());
            break;
        }

        std::map<int, bool> to_remove;
        for (std::map<int, bool>::iterator t = g1.begin(); t != g1.end(); t++)
        {
            int s = t->first;
            if (!std::isnan(record_data[s]))
            {
                to_remove[s] = 1;
            }
        }
        for (std::map<int, bool>::iterator pm = to_remove.begin(); pm != to_remove.end(); pm++)
        {
            g1.erase(pm->first);
        }
    }

    // now g1 has index of samples that are not genotyped
    for (int i = 0; i < ped.count; i++)
    {
        bool phenotyped = false;
        for (int tr = 0; tr < ped.traitNames.Length(); tr++)
        {
            if (ped[i].isPhenotyped(tr))
            {
                phenotyped = true;
                break;
            }
        }
        if (!phenotyped)
        {
            continue;
        }
        int s;
        String sample;
        if (mergedVCFID)
        {
            sample = ped[i].famid + "_" + ped[i].pid;
            s = VCFID.Integer(sample);
        } else
        {
            s = VCFID.Integer(ped[i].pid);
        }
        if (s == -1)
        {
            continue;
        }
        if (g1.find(s) != g1.end())
        {
            continue;
        }
        if (mergedVCFID)
        {
            sampleVCFIDHash.SetInteger(sample, s);
            samplePEDIDHash.SetInteger(sample, i);
        } else
        {
            sampleVCFIDHash.SetInteger(ped[i].pid, s);
            samplePEDIDHash.SetInteger(ped[i].pid, i);
        }
        genotypedSampleVCF.Push(s);
    }
}

void FastTransform::ScreenSampleID(Pedigree &ped, bool useCovariates)
{
    printf("Friendly reminder:  if your file has a lot of individuals not genotyped, the following process might take very long ...");
    fflush(stdout);
    if (PreMeta::genoFromVCF || PreMeta::dosage)
    {
        SelectSamplesVCF(ped, useCovariates);
    } else
    {
        SelectSamplesPED(ped, useCovariates);
    }
    printf("  done.\n");
}

//Prepare is to get pPheno filled up and useful persons and families clarified.
void FastTransform::Prepare(Pedigree &ped, int traitNum, bool useCovariates, FILE *log, bool shortVersion)
{
    if (!shortVersion)
    {
        printf("  Matching individuals in phenotype and genotypes ...\n");
        fprintf(log, "  Matching individuals in phenotype and genotypes ...\n");
        pPheno = new IntArray[ped.familyCount];

        //if genotpye is read from vcf file, need to make sure to eliminate individuals who are not genotyped.
        for (int p = 0; p < ped.count; p++)
        {
            if (ped[p].isFounder())
            {
                foundersHash.SetInteger(ped[p].pid, p);
            }
        }

        numFounder = 0;
        analyzedFounders = 0;
        numFamily = ped.familyCount;
        if (PreMeta::genoFromVCF || PreMeta::dosage)
        {
            for (int f = 0; f < ped.familyCount; f++)
            {
                pPheno[f].Dimension(0);
                for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                {
                    if (ped[i].isFounder())
                    {
                        numFounder++;
                    }
                    int p;
                    if (mergedVCFID)
                    {
                        String sample = ped[i].famid + "_" + ped[i].pid;
                        p = sampleVCFIDHash.Integer(sample);
                    } else
                    {
                        p = sampleVCFIDHash.Integer(ped[i].pid);
                    }

                    if (ped[i].isPhenotyped(traitNum) && (!useCovariates || ped[i].isFullyControlled()) && p != -1)
                    {
                        pPheno[f].Push(i);
                        if (ped[i].isFounder())
                        {
                            analyzedFounders++;
                        }
                    } else
                    {
                        ped[i].traits[traitNum] = _NAN_;
                    }
                }
            }
        } else
        {
            totalN = ped.count;
            // If we are using covariates, we only consider individuals
            // for which all covariates have been recorded
            for (int f = 0; f < ped.familyCount; f++)
            {
                pPheno[f].Dimension(0);
                for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                {
                    if (ped[i].isFounder())
                    {
                        numFounder++;
                    }
                    int p = samplePEDIDHash.Integer(ped[i].pid);

                    if (ped[i].isPhenotyped(traitNum) && (!useCovariates || ped[i].isFullyControlled()) && p != -1)
                    {
                        pPheno[f].Push(i);
                        if (ped[i].isFounder())
                        {
                            analyzedFounders++;
                        }
                    } else
                    {
                        ped[i].traits[traitNum] = _NAN_;
                    }
                }
            }
        }
        //printf("num of founder is: %d\n",numFounder);
        // Number of families with non-null kinships and phenotypes
        families = persons = 0;
        // Count useful families
        for (int f = 0; f < ped.familyCount; f++)
        {
            if (pPheno[f].Length())
            {
                families++;
                persons += pPheno[f].Length();
            }
        }
        if (families == 0)
        {
            printf("Trait %s has no informative families\n\n",
                   (const char *) ped.traitNames[traitNum]);
            fprintf(log, "Trait %s has no informative families\n\n",
                    (const char *) ped.traitNames[traitNum]);
        }
        printf("    Found %d phenotyped AND genotyped individuals from %d families.\n", persons, families);
        fprintf(log, "    Found %d phenotyped AND genotyped individuals from %d families.\n", persons, families);
        printf("  done.\n\n");
        fprintf(log, "  done.\n\n");
    }

    fixedEffects = 1 + (useCovariates ? ped.covariateCount : 0);

    Vector trait;
    trait.Dimension(persons);
    int index = 0;
    for (int f = 0; f < ped.familyCount; f++)
    {
        for (int j = 0; j < pPheno[f].Length(); j++)
        {
            if (ped[pPheno[f][j]].isFounder())
            {
                founders.Push(ped[pPheno[f][j]].pid);
            }
            trait[index] = ped[pPheno[f][j]].traits[traitNum];
            index++;
        }
    }

    traitVar = trait.Var();
    traitMean = trait.Average();
    UY.Dimension(persons);
    UY.Zero();
    UX.Dimension(persons, fixedEffects);
    Y.Dimension(persons);
    transU_del.Dimension(0, 0);
    // UDY.Dimension(persons);
    //  UDX.Dimension(persons,fixedEffects);
    // UD.Dimension(persons,persons);
    //inv.Dimension(fixedEffects,fixedEffects);

    //Fill in X and Y
    index = 0;
    for (int f = 0; f < ped.familyCount; f++)
    {
        int ct = pPheno[f].Length();
        for (int i = 0; i < ct; i++)
        {
            Y[index] = ped[pPheno[f][i]].traits[traitNum];
            index++;
        }
    }
    // Setup matrix of fixed effects for each family
    X.Dimension(persons, fixedEffects);
    for (int i = 0; i < persons; i++)
    { // Constant for regressing grand mean
        X[i][0] = 1.0;
    }

    if (useCovariates)
    {
        int index = 0;
        for (int f = 0; f < ped.familyCount; f++)
        {
            int count = pPheno[f].Length();
            for (int i = 0; i < count; i++)
            { // User specified covariates
                for (int j = 1; j <= ped.covariateCount; j++)
                {
                    X[index][j] = ped[pPheno[f][i]].covariates[j - 1];
                }
                index++;
            }
        }
    }
}

void FastTransform::TransformEmpkinX(Matrix &covMatrix)
{
    //printf("Transforming covariance matrix ... \n");
    //fflush(stdout);
    Matrix U;
    U.Dimension(covMatrix.rows, covMatrix.cols);
    D.Dimension(covMatrix.rows);
    EigenDecompose(covMatrix, U, D);
    transU_del.Transpose(U);
    UX.Product(transU_del, X);
    //calculate transpose(U)*y and save them in the traits vector
    for (int r = 0; r < transU_del.rows; r++)
    {
        UY[r] = transU_del[r].InnerProduct(Y);
    }
    //printf("  done.\n\n");
    FinalizeProducts();
}

//this function decompose the pedigree based var-cov matrix after VCs have been estimated

int FastTransform::GetCovarianceMatrix(AutoFit &engine, Pedigree &ped, int f, int n, Matrix &omega)
{
    if (n == -1)
    {
        return -1;
    }
    //Setup autosomal kinship contribution
    Kinship kin;
    kin.Setup(*ped.families[f]);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            omega[i][j] += 2.0 * kin(ped[pPheno[f][i]], ped[pPheno[f][j]]) * engine.heritability;
        }
    }

    //Setup X kinship contribution
    if (AutoFit::fitX && FastTransform::pedKin && !FastFit::separateX)
    {
        KinshipX kinX;
        kinX.Setup(*ped.families[f]);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                omega[i][j] += 2.0 * kinX(ped[pPheno[f][i]], ped[pPheno[f][j]]) * engine.sigma_gX;
            }
        }
    }
    //TODO:Setup shared environment contribution

    //Setup the non-shared environmentl contribution:
    //only the diagnol element should add this
    for (int i = 0; i < n; i++)
    {
        omega[i][i] += engine.sigma_e2;
    }

    /*
       printf("omega before multiply traitVar is:\n");
       for(int i=0;i<n;i++)
       {
       for(int j=0;j<n;j++)
       printf("%g ",omega[i][j]);
       printf("\n");
       }

       printf("g2,gx2,e2 are: %g,%g,%g.\n",engine.heritability,engine.sigma_gX,engine.sigma_e2);
     */
    //Finally tune-up
    omega.Multiply(engine.variance);
    return 0;
}

void FastTransform::TransformPedkinX(AutoFit &engine, Pedigree &ped)
{
    //Setup Omega
    Kinship kin;
    KinshipX kinX;
    for (int f = 0; f < ped.familyCount; f++)
    {
        int count = pPheno[f].Length();
        if (count == 0)
        {
            continue;
        }
        Matrix omega;
        omega.Dimension(count, count, 0.0);
        int status = GetCovarianceMatrix(engine, ped, f, count, omega);
        if (status == -1)
        {
            continue;
        }
        //TS:
        /*
       printf("omega is:\n");
       for(int i=0;i<count;i++)
       {
       for(int j=0;j<count;j++)
       printf("%g ",omega[i][j]);
       printf("\n");
       }
         */
        Matrix U;
        U.Dimension(count, count);
        Vector W;
        W.Dimension(count);
        EigenDecompose(omega, U, W);
        if (transU_del.rows == 0)
        {
            transU_del.Transpose(U);
            D.Copy(W);
        } else
        {
            int oldDim = transU_del.rows;
            transU_del.GrowTo(oldDim, oldDim + count, 0.0);
            transU_del.GrowTo(oldDim + count, oldDim + count, 0.0);
            D.GrowTo(transU_del.cols);
            for (int i = 0; i < count; i++)
            {
                D[i + oldDim] = W[i];
                for (int j = 0; j < count; j++)
                {
                    transU_del[i + oldDim][j + oldDim] = U[j][i];
                }
            }
        }
        /*
       printf("transU_del dimension is: %d, %d.\n",transU_del.rows,transU_del.cols);
       printf("transU_del is:\n");
       for(int i=0;i<transU_del.rows;i++)
       {
       for(int j=0;j<transU_del.cols;j++)
       printf("%g ",transU_del[i][j]);
       printf("\n");
       }
         */
    }

    //Do transformations
    UX.Product(transU_del, X);
    for (int r = 0; r < transU_del.rows; r++)
    {
        UY[r] = transU_del[r].InnerProduct(Y);
    }

    FinalizeProducts();
    /*
    //TS:
    printf("transU_del is:\n");
    for(int i=0;i<persons;i++)
    {
    for(int j=0;j<persons;j++)
    printf("%g ",transU_del[i][j]);
    printf("\n");
    }
    printf("D is:\n");
    for(int i=0;i<persons;i++)
    printf("%g\t",D[i]);
    printf("\n");
    printf("UX is:\n");
    for(int i=0;i<persons;i++)
    {
    for(int j=0;j<X.cols;j++)
    printf("%g ",UX[i][j]);
    printf("\n");
    }
     */
}

void FastTransform::TransformPedkinSepX(Pedigree &ped, int traitNum, bool useCovariates, FILE *log)
{
    //printf("Transforming pedigree kinship matrix ...\n");
    //fprintf(log,"Transforming pedigree kinship matrix ...\n");
    transU_del.Dimension(persons, persons);
    transU_del.Zero();
    D.Dimension(persons);
    int index = 0;
    for (int f = 0; f < ped.familyCount; f++)
    {
        if (pPheno[f].Length() > 0)
        {
            int count = pPheno[f].Length();
            //printf("fam %d has %d members.\n",f,count);
            Vector tmppheno;
            tmppheno.Dimension(count);
            for (int i = 0; i < count; i++)
            {
                tmppheno[i] = ped[pPheno[f][i]].traits[traitNum];
            }

            // Setup matrix of fixed effects for each family
            Matrix sub_X;
            sub_X.Dimension(count, fixedEffects);
            Matrix transU_delX;

            //fill in sub_X
            for (int i = 0; i < count; i++)
            {
                // Constant for regressing grand mean
                sub_X[i][0] = 1.0;

                // User specified covariates
                if (useCovariates)
                {
                    for (int j = 1; j <= ped.covariateCount; j++)
                    {
                        sub_X[i][j] = ped[pPheno[f][i]].covariates[j - 1];
                    }
                }
            }

            //fill in X
            if (X.rows == 0)
            {
                X.Copy(sub_X);
            } else
            {
                X.StackBottom(sub_X);
            }

            KinshipX kin;
            kin.Setup(*ped.families[f]);

            Matrix pedkin;
            pedkin.Dimension(count, count);
            for (int i = 0; i < count; i++)
            {
                for (int j = i; j < count; j++)
                {
                    pedkin[i][j] = pedkin[j][i] =
                            2.0 * kin(ped[pPheno[f][i]], ped[pPheno[f][j]]);
                }
            }
            //refine kinship matrix for adopted children
            for (int i = 0; i < count; i++)
            {
                if (ped[pPheno[f][i]].zygosity == 999)
                {
                    for (int j = 0; j < count; j++)
                    {
                        if (i != j)
                        {
                            pedkin[i][j] = pedkin[j][i] = 0.0;
                        }
                    }
                }
            }
            /*
            //TS
            printf("kinX is:\n");
            for(int i=0;i<count;i++)
            {
            for(int j=0;j<count;j++)
            printf("%g\t",pedkin[i][j]);
            printf("\n");
            }
            for(int i=0;i<count;i++)
            {
            for(int j=0;j<count;j++)
            printf("%g ",pedkin[i][j]);
            printf("\n");
            }
             */
            //now do the transformation
            //printf("Now doing the cholesky decompositon ... /n");
            //S is the vector of eigenvalues
            Vector S;
            S.Dimension(count);
            Matrix U, transU_delsub;
            U.Dimension(count, count);
            EigenDecompose(pedkin, U, S);
            transU_delsub.Transpose(U);
            transU_delX.Product(transU_delsub, sub_X);
            //fill up UX, D, and Y
            for (int row = 0; row < count; row++)
            {
                D[row + index] = S[row];
                Y[row + index] = tmppheno[row];
                for (int col = 0; col < fixedEffects; col++)
                {
                    UX[row + index][col] = transU_delX[row][col];
                }
            }

            int block = index;
            //calculate transpose(U)*y and save them in the traits vector
            for (int r = 0; r < count; r++)
            {
                //fill up transU_del here
                for (int c = 0; c < count; c++)
                {
                    transU_del[r + block][c + block] = transU_delsub[r][c];
                }
                //fill up UY here
                for (int c = 0; c < count; c++)
                {
                    UY[index] += transU_delsub[r][c] * tmppheno[c];
                }
                index++;
                //if(index>=2477 && index<=2480 || index>=2490 && index<=2492)
                //    printf("%d pid is: %s\n",index,ped[pPheno[f][r]].pid.c_str());
            }
        }
    }
    FinalizeProducts();
}

//This function will decompose pedigree based kinship for use in FAST-LMM
void FastTransform::TransformPedkin(Pedigree &ped, int traitNum, bool useCovariates, FILE *log)
{
    //printf("Transforming pedigree kinship matrix ...\n");
    //fprintf(log,"Transforming pedigree kinship matrix ...\n");
    transU_del.Dimension(persons, persons);
    transU_del.Zero();
    D.Dimension(persons);
    int index = 0;
    for (int f = 0; f < ped.familyCount; f++)
    {
        if (pPheno[f].Length() > 0)
        {
            int count = pPheno[f].Length();
            //printf("fam %d has %d members.\n",f,count);
            Vector tmppheno;
            tmppheno.Dimension(count);
            for (int i = 0; i < count; i++)
            {
                tmppheno[i] = ped[pPheno[f][i]].traits[traitNum];
            }

            // Setup matrix of fixed effects for each family
            Matrix sub_X;
            sub_X.Dimension(count, fixedEffects);
            Matrix transU_delX;

            //fill in sub_X
            for (int i = 0; i < count; i++)
            {
                // Constant for regressing grand mean
                sub_X[i][0] = 1.0;

                // User specified covariates
                if (useCovariates)
                {
                    for (int j = 1; j <= ped.covariateCount; j++)
                    {
                        sub_X[i][j] = ped[pPheno[f][i]].covariates[j - 1];
                    }
                }
            }

            //fill in X
            if (X.rows == 0)
            {
                X.Copy(sub_X);
            } else
            {
                X.StackBottom(sub_X);
            }

            Kinship kin;
            kin.Setup(*ped.families[f]);

            Matrix pedkin;
            pedkin.Dimension(count, count);
            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    pedkin[i][j] = 2.0 * kin(ped[pPheno[f][i]], ped[pPheno[f][j]]);
                }
            }
            //refine kinship matrix for adopted children
            for (int i = 0; i < count; i++)
            {
                if (ped[pPheno[f][i]].zygosity == 999)
                {
                    for (int j = 0; j < count; j++)
                    {
                        if (i != j)
                        {
                            pedkin[i][j] = pedkin[j][i] = 0.0;
                        }
                    }
                }
            }
            //now do the transformation
            //printf("Now doing the cholesky decompositon ... /n");
            //S is the vector of eigenvalues
            Vector S;
            S.Dimension(count);
            Matrix U, transU_delsub;
            U.Dimension(count, count);
            EigenDecompose(pedkin, U, S);
            transU_delsub.Transpose(U);
            transU_delX.Product(transU_delsub, sub_X);

            //fill up UX, D, and Y
            for (int row = 0; row < count; row++)
            {
                D[row + index] = S[row];
                Y[row + index] = tmppheno[row];
                for (int col = 0; col < fixedEffects; col++)
                {
                    UX[row + index][col] = transU_delX[row][col];
                }
            }

            int block = index;
            //calculate transpose(U)*y and save them in the traits vector
            for (int r = 0; r < count; r++)
            {
                //fill up transU_del here
                for (int c = 0; c < count; c++)
                {
                    transU_del[r + block][c + block] = transU_delsub[r][c];
                }
                //fill up UY here
                for (int c = 0; c < count; c++)
                {
                    UY[index] += transU_delsub[r][c] * tmppheno[c];
                }
                index++;
                //if(index>=2477 && index<=2480 || index>=2490 && index<=2492)
                //    printf("%d pid is: %s\n",index,ped[pPheno[f][r]].pid.c_str());
            }
        }
    }
    //printf(" done.\n");
    //fprintf(log," done.\n");
    FinalizeProducts();
}

//This function will decompose kinship for use in FAST-LMM
void FastTransform::TransformEmpkin(Pedigree &ped, int traitNum, bool useCovariates, KinshipEmp &kin_emp, Matrix &input,
                                    FILE *log)
{
    transU_del.Dimension(persons, persons);
    D.Dimension(persons);
    Matrix U;
    U.Dimension(persons, persons);
    Matrix allPairs;
    allPairs.Dimension(persons, persons);

    //if reading in Empirical matrix from outside
    if (readInEmp != "")
    {
        printf("\n    Matching IDs in kinship matrix ... ");
        fprintf(log, "\n    Matching IDs in kinship matrix ... ");
        fflush(stdout);
        IntArray selectedSample;
        //Checking if all individuals are included in read-in kinship.
        for (int f = 0; f < ped.familyCount; f++)
        {
            for (int j = 0; j < pPheno[f].Length(); j++)
            {
                int p;
                if (mergedVCFID)
                {
                    p = kin_emp.IDFromEmp.Integer(ped[pPheno[f][j]].famid + "_" + ped[pPheno[f][j]].pid);
                } else
                {
                    p = kin_emp.IDFromEmp.Integer(ped[pPheno[f][j]].pid);
                }

                if (p == -1)
                {
                    error("Person %s in ped file is not included in read-in kinship.\n", ped[pPheno[f][j]].pid.c_str());
                    fprintf(log, "Person %s in ped file is not included in read-in kinship.\n",
                            ped[pPheno[f][j]].pid.c_str());
                } else
                {
                    selectedSample.Push(p);
                }
            }
        }

        //subsetting
        for (int r = 0; r < selectedSample.Length(); r++)
        {
            for (int c = 0; c < selectedSample.Length(); c++)
            {
                allPairs[r][c] = input[selectedSample[r]][selectedSample[c]];
            }
        }

        printf("    done.\n");
        fprintf(log, "    done.\n");
    } else
    {
        //printf("Calculating empirical kinship matrix ...\n");
        printf("\n    Matching IDs in kinship matrix ... ");
        fprintf(log, "\n    Matching IDs in kinship matrix ... ");
        fflush(stdout);
        SubSet(input, ped, allPairs);
        printf("done.\n");
        fprintf(log, "done.\n");
    }

    printf("    Decomposing empirical kinship matrix ... ");
    fprintf(log, "    Decomposing empirical kinship matrix ... ");
    fflush(stdout);
    EigenDecompose(allPairs, U, D);
    printf("done.\n");
    fprintf(log, "done.\n");

    transU_del.Transpose(U);
    UX.Product(transU_del, X);
    //calculate transpose(U)*y and save them in the traits vector
    for (int r = 0; r < transU_del.rows; r++)
    {
        UY[r] = transU_del[r].InnerProduct(Y);
    }
    FinalizeProducts();
}

void FastTransform::SubSet(Matrix &allPairs, Pedigree &ped, Matrix &tmp)
{
    IntArray sampleOrdered;
    if (PreMeta::genoFromPed)
    {
        for (int f = 0; f < ped.familyCount; f++)
        {
            for (int i = 0; i < pPheno[f].Length(); i++)
            {
                for (int idx = 0; idx < genotypedSamplePED.Length(); idx++)
                {
                    if (pPheno[f][i] == genotypedSamplePED[idx])
                    {
                        sampleOrdered.Push(idx);
                        break;
                    }
                }
            }
        }
    }

    if (PreMeta::genoFromVCF || PreMeta::dosage)
    {
        for (int f = 0; f < ped.familyCount; f++)
        {
            for (int j = 0; j < pPheno[f].Length(); j++)
            {
                int p;
                if (mergedVCFID)
                {
                    String sample = ped[pPheno[f][j]].famid + "_" + ped[pPheno[f][j]].pid;
                    p = sampleVCFIDHash.Integer(sample);
                } else
                {
                    p = sampleVCFIDHash.Integer(ped[pPheno[f][j]].pid);
                }
                for (int idx = 0; idx < genotypedSampleVCF.Length(); idx++)
                {
                    if (p == genotypedSampleVCF[idx])
                    {
                        sampleOrdered.Push(idx);
                        break;
                    }
                }
            }
        }
    }

    //printf("sampleOrdered length: %d; persons: %d\n",sampleOrdered.Length(),persons);
    for (int r = 0; r < sampleOrdered.Length(); r++)
    {
        for (int c = 0; c < sampleOrdered.Length(); c++)
        {
            tmp[r][c] = allPairs[sampleOrdered[r]][sampleOrdered[c]];
        }
    }
}

int FastTransform::EigenDecompose(Matrix &matrix_in, Matrix &U, Vector &D)
{
    if (matrix_in.rows == 0 || matrix_in.cols == 0)
    {
        return -1;
    }
    if (U.rows == 0 || U.cols == 0)
    {
        return -1;
    }
    if (D.Length() == 0)
    {
        return -1;
    }
    if (matrix_in.rows != matrix_in.cols)
    {
        error("Your kinship matrix is not symmetric!\n");
    }

    int n = matrix_in.rows;
    Eigen::MatrixXf todo(n, n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            todo(i, j) = matrix_in[i][j];
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(todo);

    if (solver.info() == Eigen::Success)
    { //copy eigen vectors and eigenvalues
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                U[i][j] = solver.eigenvectors()(i, j);
            }
        }
        for (int i = 0; i < n; i++)
        {
            D[i] = solver.eigenvalues()[i];
            if (D[i] < 0)
            {
                for (int j = 0; j < n; j++)
                {
                    U[j][i] *= -1.0;
                }
                D[i] *= -1.0;
            }
        }
        return 0;
    }
    return -1;
}

bool FastTransform::FinalizeProducts()
{
    if (transU_del.rows != persons || transU_del.cols != persons)
    {
        return false;
    }
    transU.resize(persons, persons);
    //copy transU_del to transU and delete transU_del
    for (int i = 0; i < transU_del.rows; i++)
    {
        for (int j = 0; j < transU_del.cols; j++)
        {
            transU(i, j) = transU_del[i][j];
        }
    }
    transU_del.Dimension(0, 0);
    return true;
}

