#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"
#include "driftkickrad.c"	/* drift6.c, bndthinkickrad.c */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656


#define SQR(X) ((X)*(X))

struct elem 
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    double BendingAngle;
    double EntranceAngle;
    double ExitAngle;
    double Energy;
    /* Optional fields */
    int FringeBendEntrance;
    int FringeBendExit;
    double FringeInt1;
    double FringeInt2;
    double FullGap;
    int FringeQuadEntrance;
    int FringeQuadExit;
    double *fringeIntM0;
    double *fringeIntP0;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
};

void BndMPoleSymplectic4RadPass(double *r, double le, double irho, double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, 	double exit_angle,
        int FringeBendEntrance, int FringeBendExit,
        double fint1, double fint2, double gap,
        int FringeQuadEntrance, int FringeQuadExit,
        double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
        double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */       
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double E0, int num_particles)
        
{	int c,m;
    double *r6;
    double SL, L1, L2, K1, K2;
    bool useT1, useT2, useR1, useR2, useFringe1, useFringe2;
    bool useLinFrEleEntrance = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadEntrance==2);
    bool useLinFrEleExit = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadExit==2);
    SL = le/num_int_steps;
    L1 = SL*DRIFT1;
    L2 = SL*DRIFT2;
    K1 = SL*KICK1;
    K2 = SL*KICK2;
    if (T1==NULL)
        useT1=false;
    else
        useT1=true;
    if (T2==NULL)
        useT2=false;
    else
        useT2=true;
    if (R1==NULL)
        useR1=false;
    else
        useR1=true;
    if (R2==NULL)
        useR2=false;
    else
        useR2=true;
    
    /* calculate entrance fringe field if fint, gap and FringeBendEntrance are not 0 */
    if( fint1==0 || gap==0 || FringeBendEntrance==0)
        useFringe1 = false;
    else
        useFringe1=true;
    /* calculate exit fringe field if fint, gap and FringeBendExit are not 0 */
    if( fint2==0 || gap==0 || FringeBendExit==0)
        useFringe2 = false;
    else
        useFringe2=true;
    
    for(c = 0;c<num_particles;c++)	/* Loop over particles */
    {   
        r6 = r+c*6;
        if(!atIsNaN(r6[0]))
        {
            /*  misalignment at entrance  */
            if (useT1)
                ATaddvv(r6,T1);
            if (useR1)
                ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* edge focus */
             if(useFringe1)
                if (FringeBendEntrance==1)
                    edge_fringe(r6, irho, entrance_angle,fint1,gap);
                else if (FringeBendEntrance==2)
                    edge_fringe_Version2(r6, irho, entrance_angle,fint1,gap);
                else
                    edge_fringe_Version3Entrance(r6, irho, entrance_angle,fint1,gap);
            else
                edge(r6, irho, entrance_angle);
            /* quadrupole gradient fringe */
            if (FringeQuadEntrance && B[1]!=0)
                if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantEntrance(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassP(r6,B[1]);
            /* integrator  */
            for(m=0; m < num_int_steps; m++) /* Loop over slices */
            {
                drift6(r6,L1);
                bndthinkickrad(r6, A, B, K1, irho, E0, max_order);
                drift6(r6,L2);
                bndthinkickrad(r6, A, B, K2, irho, E0, max_order);
                drift6(r6,L2);
                bndthinkickrad(r6, A, B,  K1, irho, E0, max_order);
                drift6(r6,L1);
            }
            /* quadrupole gradient fringe */
            if (FringeQuadExit && B[1]!=0)
                if (useLinFrEleExit) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantExit(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassN(r6,B[1]);
             /* edge focus */
           if(useFringe2)
                if (FringeBendExit==1)
                    edge_fringe(r6, irho, entrance_angle,fint1,gap);
                else if (FringeBendExit==2)
                    edge_fringe_Version2(r6, irho, entrance_angle,fint1,gap);
                else
                    edge_fringe_Version3Exit(r6, irho, entrance_angle,fint1,gap);
            else
                edge(r6, irho, exit_angle);
            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* Misalignment at exit */
            if (useR2)
                ATmultmv(r6,R2);
            if (useT2)
                ATaddvv(r6,T2);
        } 
    }
}



#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    double irho;
    if (!Elem) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap,
                FringeInt1, FringeInt2, Energy;
        int MaxOrder, NumIntSteps,  FringeBendEntrance, FringeBendExit,
                FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",1); check_error();
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",1); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        Elem->BendingAngle=BendingAngle;
        Elem->EntranceAngle=EntranceAngle;
        Elem->ExitAngle=ExitAngle;
        Elem->Energy=Energy;
        /*optional fields*/
        Elem->FullGap=FullGap;
        Elem->FringeInt1=FringeInt1;
        Elem->FringeInt2=FringeInt2;
        Elem->FringeBendEntrance=FringeBendEntrance;
        Elem->FringeBendExit=FringeBendExit;
        Elem->FringeQuadEntrance=FringeQuadEntrance;
        Elem->FringeQuadExit=FringeQuadExit;
        Elem->fringeIntM0=fringeIntM0;
        Elem->fringeIntP0=fringeIntP0;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
    }
    irho = Elem->BendingAngle/Elem->Length;
    BndMPoleSymplectic4RadPass(r_in,Elem->Length,irho,Elem->PolynomA,Elem->PolynomB,
            Elem->MaxOrder,Elem->NumIntSteps,Elem->EntranceAngle,Elem->ExitAngle,
            Elem->FringeBendEntrance,Elem->FringeBendExit,
            Elem->FringeInt1,Elem->FringeInt2,Elem->FullGap,
            Elem->FringeQuadEntrance,Elem->FringeQuadExit,
            Elem->fringeIntM0,Elem->fringeIntP0,Elem->T1,Elem->T2,
            Elem->R1,Elem->R2,Elem->RApertures,Elem->EApertures,Elem->Energy,num_particles);
    return Elem;
}

MODULE_DEF(BndMPoleSymplectic4Pass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, 
                FringeInt1, FringeInt2, Energy;
        int MaxOrder, NumIntSteps, FringeBendEntrance, FringeBendExit,
                FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        double irho;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        irho = BendingAngle/Length;
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        BndMPoleSymplectic4RadPass(r_in, Length, irho, PolynomA, PolynomB,
                MaxOrder,NumIntSteps,EntranceAngle,ExitAngle,
                FringeBendEntrance,FringeBendExit,FringeInt1,FringeInt2,
                FullGap,FringeQuadEntrance,FringeQuadExit,fringeIntM0,fringeIntP0,
                T1,T2,R1,R2,RApertures,EApertures,Energy,num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
        mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
        mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],5,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],6,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],7,mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0],8,mxCreateString("Energy"));
        
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(15,1);
            mxSetCell(plhs[1],0,mxCreateString("FullGap"));
            mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
            mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
            mxSetCell(plhs[1],3,mxCreateString("FringeBendEntrance"));
            mxSetCell(plhs[1],4,mxCreateString("FringeBendExit"));
            mxSetCell(plhs[1],5,mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1],6,mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1],7,mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1],8,mxCreateString("fringeIntP0"));
            mxSetCell(plhs[1],9,mxCreateString("T1"));
            mxSetCell(plhs[1],10,mxCreateString("T2"));
            mxSetCell(plhs[1],11,mxCreateString("R1"));
            mxSetCell(plhs[1],12,mxCreateString("R2"));
            mxSetCell(plhs[1],13,mxCreateString("RApertures"));
            mxSetCell(plhs[1],14,mxCreateString("EApertures"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */


