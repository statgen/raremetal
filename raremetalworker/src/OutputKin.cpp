#include "OutputKin.h"
#include "TransformResiduals.h"
#include "PreMeta.h"

bool OutputKin::outputKinX = false;
bool OutputKin::outputKin = false;
bool OutputKin::kinOnly = false;

/*
void OutputKin::Write(Family & f,int before,int N)
{  
   FILE * file = NULL;
   String filename;
   filename.printf("Pedigree-Kinship-%s",(const char *) PreMeta::outputFile);
   file = fopen(filename, "at");
   Kinship kin;
   kin.Setup(f);

   int famNum = f.count;
   int after = N-famNum-before;

   if (file != NULL)
      for(int i=0;i<famNum;i++)
      {
	 for(int j=0;j<before;j++)
	    fprintf(file,"0\t");
	 for(int k=0;k<famNum;k++)
	    fprintf(file,"%g\t",kin.allPairs[i][k]);
	 for(int l=0;l<after;l++)
	    fprintf(file,"0\t");
	 fprintf(file,"\n");
      }

   if (file != 0)
   {
      fclose(file);
   }
}

void OutputKin::WriteX(Family & f,int before,int N)
{  
   FILE * file = NULL;
   String filename;
   filename.printf("KinshipX-%s",(const char *) PreMeta::outputFile);
   file = fopen(filename, "wt");
   KinshipX kin;
   kin.Setup(f);
   int famNum = f.count;
   int after = N-famNum-before;
   if (file != NULL)
      for(int i=0;i<famNum;i++)
      {
	 for(int j=0;j<before;j++)
	    fprintf(file,"0\t");
	 for(int k=0;k<famNum;k++)
	    fprintf(file,"%g\t",kin.allPairs[i][k]);
	 for(int l=0;l<after;l++)
	    fprintf(file,"0\t");
	 fprintf(file,"\n");
      }
   if (file != 0)
   {
      fclose(file);
   }
}

void OutputKin::Output(Pedigree & ped)
{
   if(outputKin && !FastTransform::empKin)
   {
      //output pedigree kinship matrix
      int N = ped.count;
      int counter = 0;
      for(int i=0; i<ped.familyCount; ++i)
      {
	 Family *f = ped.families[i];
	 int famNum = f->count;
	 Write(*f,counter,N);
	 counter += famNum;
      }
   }

   if(outputKinX && !FastTransform::empKin)
   {
      //output pedigree kinshipX matrix
      int N = ped.count;
      int counter = 0;
      for(int i=0; i<ped.familyCount; ++i)
      {
	 Family *f = ped.families[i];
	 int famNum = f->count;
	 WriteX(*f,counter,N);
	 counter += famNum;
      }
   }
}
*/
